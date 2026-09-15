#include <gmp.h>
#include <cmath>
#include <DMCrand.h>
#include <bit_array.h>
#include <map>
#include <algorithm>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <atomic>
#include <future>
#include <ctime>
#include <limits>
#include <sys/resource.h>

using namespace std;

struct mpf_pair {
    mpf_t* Pr;
    mpf_t z;

    mpf_pair(mpf_t *pr, const mpf_t& zz) : Pr(pr) {
        mpf_init(z);
        mpf_set(z, zz);
    }

    ~mpf_pair() {
        mpf_clear(z);
    }
};

mpf_pair prob_CG(int*& dnf, long*& dnf_ind, int& m, int& n) {
    mpf_t* Pr = new mpf_t[m];
    mpf_t z;
    mpf_init(z);
    mpf_set_d(z,0.0);

    for (int j=0; j<m; ++j) {

        mpf_init(Pr[j]);
        int size = dnf[dnf_ind[j]];

        mpf_set_d(Pr[j], 1/pow(2.0,size));
        mpf_add(z,z,Pr[j]);
    }

    return {Pr, z};
}

// Succinct sampling from discrete distributions, Bringmann and Larsen
pair<int*, int> sample_info(mpf_t*& Pr, int& m, mpf_t& z) {
    int* A = new int[2*m]();

    int loc = 0;
    int count = 0;
    int counts = 0;

    mpf_t temp;
    mpf_init(temp);

    for (int i=0; i<m; ++i) {
        mpf_mul_ui(temp, Pr[i], m);
        mpf_div(temp, temp, z);
        mpf_floor(temp,temp);
        counts = mpf_get_si(temp);

        count = counts + 1;

        for (int k=0; k < count; ++k) {
            A[loc+k] = i;
        }

        loc += count;
    }

    return {A, loc};
}

// Succinct sampling from discrete distributions, Bringmann and Larsen
int select_clause(int& u, int*& A, mpf_t*& Pr, mpf_t& z, int& nz, int& m, mt19937& gen, 
        uniform_real_distribution<double>& dis, uniform_int_distribution<int>& idis) {

    u = -1;
    while (u == -1) {
        u = uniform_integer(nz, gen);

        // rejection algorithm
        if (u == 0 || A[u-1] != A[u]) {
            static mpf_t temp_int;
            static mpf_t fractional_mpf;
            static mpf_t prob;
            static int init = 0;

            if (!init) {
                init = 1;   
                mpf_init(temp_int);
                mpf_init(fractional_mpf);
                mpf_init(prob);
            }

            mpf_set_d(prob, static_cast<double>(m));
            mpf_mul(prob, prob, Pr[A[u]]);
            mpf_div(prob, prob, z);

            mpf_floor(temp_int, prob);
            mpf_sub(fractional_mpf, prob, temp_int);

            double fractional = mpf_get_d(fractional_mpf);

            double rej = generate_random_real(gen,dis);
            if (rej >= (fractional)) {
                u = -1;
            }
        }
    }

    return A[u];
}

void lazyassign(int*& dnf, long& C_ind, int& n, int& bytesize, gmt::byte*& clause_vars,
        gmt::byte*& set_vars, bit_array_write<gmt::byte, int>& bitWrite, bit_array_read<gmt::byte, int>& bitRead) {

    for (int j=1; j<dnf[C_ind]+1; ++j) {
        int var = dnf[C_ind+j];
        int varI = abs(var)-1;
        bitWrite(clause_vars, varI);
        if (var > 0) {
            bitWrite(set_vars, varI);
        } else {
            bitWrite.clear(set_vars, varI);
        }
    }
}

bool lazySAT(int*& dnf, long& C_ind, gmt::byte* const& set_vars, gmt::byte*& clause_vars,
        bit_array_read<gmt::byte, int>& bitRead, bit_array_write<gmt::byte, int>& bitWrite,
        mt19937_64& gen, uniform_int_distribution<unsigned long>& dist, unsigned long& buff, int& buff_length) {
    const int size = dnf[C_ind] + 1;
    int var;
    int varI;
    int sign;

    for (int j=1; j<size; ++j) {
        var = dnf[C_ind+j];
        varI = abs(var)-1;

        if (!bitRead(clause_vars, varI)) {
            bitWrite(clause_vars, varI);

            sign = true_or_false(gen, dist, buff, buff_length, sign);
            if (sign > 0) {
                bitWrite(set_vars, varI);
                if (var < 0) {
                    return false;
                }
            } else {
                bitWrite.clear(set_vars, varI);
                if (var > 0) {
                    return false;
                }
            }
        } else if ((var > 0) ^ (bitRead(set_vars, varI))) {
            return false;
        } 
    }

    return true;
}

// read a .DNF file. first line specifies "p dnf n m" where n,m are number of variables and clauses
// each clause is on a new line and ends with a 0
tuple<int**, int, int, unsigned long> readDNF(const string& filename) {
    ifstream file(filename);
    string line;

    int count = 0;
    int n;
    int m;

    if (getline(file,line)) {
        istringstream iss(line);
        string info;
        while (iss >> info) {
            count += 1;
            if (count == 3) {
                n = stoi(info);
            }
            if (count == 4) {
                m = stoi(info);
                break;
            }
        }
    }

    int** T = new int*[m];
    int* clause = new int[200];
    unsigned long w_tot = 0;

    int clausecount = 0;
    while (getline(file, line)) {
        istringstream iss(line);
        int value;
        count = 0;

        while (iss >> value) {

            if (value == 0) {

                int* clause_compress = new int[count+1];
                clause_compress[0] = count;
                memcpy(clause_compress+1,clause,count*sizeof(int));

                T[clausecount] = clause_compress;
                clausecount += 1;
                w_tot += count;
                break;
            } else {
                clause[count] = value;
                count += 1;
            }
        }
    }

    delete[] clause;

    file.close();
    return make_tuple(T, n, m, w_tot);
}

void timeout(atomic<bool>& running,
             int cutoff_seconds,
             condition_variable& cv,
             mutex& mtx)
{
    auto deadline = chrono::steady_clock::now()
                    + chrono::seconds(cutoff_seconds);

    unique_lock<mutex> lock(mtx);
    cv.wait_until(lock, deadline,
                  [&] { return !running.load(); });

    if (running.exchange(false)) {
        cv.notify_all();
    }
}

double cpuTime(void) {
    struct rusage ru;
    getrusage(RUSAGE_SELF,&ru);
    return (double)ru.ru_utime.tv_sec + (double)ru.ru_utime.tv_usec / 1000000.0;
}