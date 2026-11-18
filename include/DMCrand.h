#ifndef DMCrand_H
#define DMCrand_H

#include <random>

using namespace std;

double generate_random_real(mt19937& gen, uniform_real_distribution<double>& dis);

// int uniform_integer(uniform_int_distribution<int>& dis, mt19937& gen);

int uniform_integer(int m, mt19937& gen);

int weighted_torf(mt19937& gen, uniform_real_distribution<double>& dis, double weight);

int true_or_false(mt19937_64& gen, uniform_int_distribution<unsigned long>& dist, unsigned long& buff, int& buff_length, int& sign);

#endif