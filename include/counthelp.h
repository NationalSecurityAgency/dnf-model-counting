#ifndef DMChelp_H
#define DMChelp_H

#include <gmp.h>
#include <cmath>
#include <DMCrand.h>
#include <bit_array.h>
#include <map>
#include <atomic>
#include <future>

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

mpf_pair prob_CG(int*& dnf, long*& dnf_ind, int& m, int& n);

pair<int*, int> sample_info(mpf_t*& Pr, int& m, mpf_t& z);

int select_clause(int& u,int*& A, mpf_t*& Pr, mpf_t& z, int& nz, int& m, mt19937& gen, uniform_real_distribution<double>& dis, uniform_int_distribution<int>& idis);

void lazyassign(int*& dnf, long& C_ind, int& n, int& bytesize, gmt::byte*& clause_vars,
        gmt::byte*& set_vars, bit_array_write<gmt::byte, int>& bitWrite, bit_array_read<gmt::byte, int>& bitRead);

bool lazySAT(int*& dnf, long& C_ind, gmt::byte* const& set_vars, gmt::byte*& clause_vars,
        bit_array_read<gmt::byte, int>& bitRead, bit_array_write<gmt::byte, int>& bitWrite,
        mt19937_64& gen, uniform_int_distribution<unsigned long>& dist, unsigned long& buff, int& buff_length);

tuple<int**, int, int, unsigned long> readDNF(const string& filename);

void timeout(atomic<bool>& running, int timeout_seconds, condition_variable& cv, mutex& mtx);

double cpuTime(void);

#endif