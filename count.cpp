#include <DMCrand.h>
#include <random>
#include <fstream>
#include <cstring>
#include <iostream>
#include <bit_array.h>
#include <DMChelp.h>
#include <map>
#include <sys/time.h>
#include <gmp.h>
#include <chrono>
#include <thread>
#include <atomic>
#include <future>
#include <algorithm>
#include <ctime>
#include <cmath>

using namespace std;

void shuffle(int**& DNF, long*& dnf_ind, int*& dnf, int*& dnfMax, int*& set_vars, int& current_var ,int& level, int& m, 
    unsigned long& form_size, unsigned long& current_loc, 
    unsigned long& w_tot, int& w_min, double& beta, double& rand_clause, int& min_loc) {
    int swap;
    rand_clause = drand48();
    int* k = 0;
    int clause_size;
    unsigned long clause_loc;
    
    int var = 0;
    int sign = 1;

    // decide whether to swap with a random clause or with the least remaining ordered clause
    if (rand_clause < 1.0-(1.0-beta)*min(1.0,w_tot/((m-level)*static_cast<double>(w_min)))) {
        while (true) {
            swap = rand() % (m - min_loc) + min_loc; // random remaining clause
            k = DNF[swap];
            if (k!=0) {
                break;
            }
        }

        clause_loc = form_size - current_loc;
        dnf_ind[level] = clause_loc;
        DNF[swap] = 0;

        clause_size = k[0];
        dnf[clause_loc] = clause_size;
        w_tot -= clause_size;

        // remap variables for tracking max var up to this clause
        for (int i=1; i<clause_size+1; ++i) {
            var = abs(k[i]);
            if (k[i] > 0) {
                sign = 1;
            } else {
                sign = -1;
            }
            if (set_vars[var] == 0) {
                set_vars[var] = current_var;
                dnf[clause_loc + i] = sign * current_var;
                ++current_var;
            } else {
                dnf[clause_loc + i] = sign * set_vars[var];
            }
        }

        // sort clause vars and store max var
        sort(dnf + clause_loc + 1, dnf + clause_loc + clause_size + 1, [](int a, int b) {return abs(a) < abs(b);});
        dnfMax[level] = current_var-1;
        current_loc -= k[0]+1;

        // calculate least remaining clause
        if (swap == min_loc) {
            while (true) {
                min_loc += 1;
                k = DNF[min_loc];
                if (k!=0) {
                    break;
                }
            }
            w_min = DNF[min_loc][0];
        }
    } else {
        k = DNF[min_loc]; // least remaining clause

        clause_loc = form_size - current_loc;
        dnf_ind[level] = clause_loc;
        DNF[min_loc] = 0;

        clause_size = k[0];
        dnf[clause_loc] = clause_size;
        w_tot -= k[0];

        // remap variables for tracking max var up to this clause
        for (int i=1; i<clause_size+1; ++i) {
            var = abs(k[i]);
            if (k[i] > 0) {
                sign = 1;
            } else {
                sign = -1;
            }
            if (set_vars[var] == 0) {
                set_vars[var] = current_var;
                dnf[clause_loc + i] = sign * current_var;
                ++current_var;
            } else {
                dnf[clause_loc + i] = sign * set_vars[var];
            }
        }
        // sort clause vars and store max var
        sort(dnf + clause_loc + 1, dnf + clause_loc + clause_size + 1, [](int a, int b) {return abs(a) < abs(b);});
        current_loc -= clause_size+1;
        dnfMax[level] = current_var-1;

        // calculate least remaining clause
        while(true) {
            min_loc += 1;
            k = DNF[min_loc];
            if (k!=0) {
                break;
            }
        }
        w_min = DNF[min_loc][0];
    }
}

double ouralg(int*& dnf, long*& dnf_ind, int*& dnfMax, int& maxvar, int*& A, mpf_t*& Pr, int& nz,
    mpf_t& z, int& nvars, int& m, double& eps, double& delta, bit_array_read<gmt::byte, int>& bitRead,
    bit_array_write<gmt::byte, int>& bitWrite, mt19937& gen, mt19937_64& gen64, uniform_real_distribution<double>& dis, uniform_int_distribution<int>& idis, int Thresh) {

    int Y = 0;
    long count = 0;

    int bytesize = int(ceil(maxvar/8.0));
    // set_vars stores the affirmation or negation of the variable as a 1 or 0 bit
    gmt::byte* set_vars = new gmt::byte[bytesize];
    // clause_vars stores which variables have been assigned
    gmt::byte* clause_vars = new gmt::byte[bytesize]();

    uniform_int_distribution<unsigned long> dist(numeric_limits<unsigned long>::min(), numeric_limits<unsigned long>::max());
    int index;
    long test_ind;
    long C_ind;
    bool clause_true;
    unsigned long buff = 0;
    int buff_length = 0;

    // Trial loop
    while (Y < Thresh) {
        // Select a clause with probability proportional to its weight
        index = select_clause(index, A, Pr, z, nz, m, gen, dis, idis);
        C_ind = dnf_ind[index];

        int max_var = -1;

        // assign only variables from the selected clause
        lazyassign(dnf, C_ind, nvars, bytesize, clause_vars, set_vars, bitWrite, bitRead);

        double t1 = 1;
        bool abort = false;
        double u1 = generate_random_real(gen,dis);

        // step loop
        for (int level=0; level<m; ++level) {

            if (level == index) {
                continue;
            }

            test_ind = dnf_ind[level];
            clause_true = lazySAT(dnf, test_ind, set_vars, clause_vars, bitRead, bitWrite, gen64, dist, buff, buff_length);

            if (clause_true == true) {
                t1 += 1;

                if (t1 >= 1/u1) {
                    abort = true;
                    max_var = dnfMax[level];
                    break;
                }
            }
        }

        if (abort == false) {
            Y += 1;
            max_var = maxvar;
        }
        for (int i=dnf[C_ind]; i>0; --i) {
            if (abs(dnf[C_ind+i])-1 <= max_var) {
                break;
            }
            bitWrite.clear_byte(clause_vars,abs(dnf[C_ind+i])-1);
        }
        memset(clause_vars, 0, ceil((max_var)/8.0));
        count += 1;
    }

    delete[] set_vars;
    delete[] clause_vars;

    return 1.0*Y/count;
}

int main() {
    int exp = 0;
    int nvars = 4;
    double eps = .05;
    double delta = .05;

    mpf_set_default_prec(128);

    bit_array_read<gmt::byte,int> bitRead;
    bit_array_write<gmt::byte,int> bitWrite;

    srand(time(0));
    random_device rd;
    mt19937 gen(rd());
    mt19937_64 gen64(rd());

    tuple<int**, int, int, unsigned long> DNFINFO = readDNF("test.dnf");

    int** DNF = get<0>(DNFINFO);
    int n = get<1>(DNFINFO);
    int m = get<2>(DNFINFO); 
    unsigned long w_tot = get<3>(DNFINFO);

    uniform_real_distribution<double> dis(0.0,1.0);

    double beta = 0.01;
    unsigned long form_size = w_tot + m;

    timeval start;
    timeval end;
    gettimeofday(&start,0);

    double rand_clause;

    // create the ordered partition
    sort(DNF, DNF+m, [](int* a, int* b) {return a[0] < b[0];});

    int w_min = DNF[0][0];

    long* dnf_ind = new long[m];
    int* dnfMax = new int[m];
    int* dnf = new int[form_size];

    int min_loc = 0;

    unsigned long current_loc = form_size;
    // during the shuffling, set_vars will remap unique vars as encountered to keep track of the maximum variable encountered up to a given clause.
    int* set_vars = new int[n+1]();
    int current_var = 1;

    unsigned long randN = 0;

    for (int j=0; j<m-1; ++j) {
        shuffle(DNF, dnf_ind, dnf, dnfMax, set_vars, current_var, j, m, form_size, current_loc, w_tot, w_min, beta, rand_clause, min_loc);
    }

    unsigned long clause_loc = form_size - current_loc;
    dnf_ind[m-1] = clause_loc;
    int* k = DNF[min_loc];
    int clause_size = k[0];
    dnf[clause_loc] = clause_size;
    int var = 0;
    int sign = 1;

    for (int i=1; i<clause_size+1; ++i) {
        var = abs(k[i]);
        if (k[i] > 0) {
            sign = 1;
        } else {
            sign = -1;
        }
        if (set_vars[var] == 0) {
            set_vars[var] = current_var;
            dnf[clause_loc + i] = sign * current_var;
            ++current_var;
        } else {
            dnf[clause_loc + i] = sign * set_vars[var];
        }
    }
    dnfMax[m-1] = --current_var;

    // Calculate the threshold T for the number of successful trials needed
    int T = ceil( log(2/delta) / (log(1+eps) - eps/(1+eps)) );

    double upper = std::exp(eps/(1+eps)) / (1+eps);
    double lower = std::exp(-eps/(1-eps)) / (1-eps);
    while (pow(upper,T) + pow(lower,T) <= delta) {
        T--;
    }
    T++;

    // Compute the clause selection distribution
    mpf_pair prob_table1 = prob_CG(dnf, dnf_ind, m, n);
    mpf_t* Pr1 = prob_table1.Pr;
    mpf_t z_mpf1;
    mpf_init(z_mpf1);
    mpf_set(z_mpf1, prob_table1.z);
    double z1 = mpf_get_d(z_mpf1);

    // Form the data structure for efficient sampling
    pair<int*, int> SAMPLING1 = sample_info(Pr1, m, z_mpf1);
    int* A1 = SAMPLING1.first;
    int nz1 = SAMPLING1.second;
    uniform_int_distribution<int> idis(0,nz1-1);


    double overlap = ouralg(dnf, dnf_ind, dnfMax, current_var, A1, Pr1, nz1, z_mpf1, n, m, eps, delta, bitRead, bitWrite, gen, gen64, dis, idis, T);

    double approx = overlap * z1;

    cout << approx << endl;

    return 0;
}

// g++ -lgmp -I include/ -O3 -o count count.cpp rand.cpp counthelp.cpp

