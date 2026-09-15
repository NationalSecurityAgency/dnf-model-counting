#include <random>

using namespace std;

double generate_random_real(mt19937& gen, uniform_real_distribution<double>& dis) {
    return dis(gen);
}

int uniform_integer(int m, mt19937& gen) {
    return rand() % m;
}

// assign variables true or false with probability determined by its weight
int weighted_torf(mt19937& gen, uniform_real_distribution<double>& dis, double weight) {
    double u = dis(gen);

    if (u <= weight) {
        return 1;
    } else {
        return -1;
    }
}

// assign variables true or false using a 64 bit number
int true_or_false(mt19937_64& gen, uniform_int_distribution<unsigned long>& dist, unsigned long& buff, int& buff_length, int& sign) {

    if (buff_length == 0) {
        buff = dist(gen);
        buff_length = 64;
    }

    buff_length--;

    sign = buff & 1;
    buff = buff >> 1;

    return sign;
}