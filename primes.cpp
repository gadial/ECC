#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "primes.h"

#include <stdio.h>
#define MILLER_RABIN_REPEATS 5

RandomNumberGenerator::RandomNumberGenerator(unsigned long seed){
    gmp_randinit_mt(state);
    if (seed == 0)
        seed = time(NULL);
    gmp_randseed_ui(state,seed);
}

mpz_class RandomNumberGenerator::rand(mpz_class max){
    mpz_t temp;
    mpz_init(temp);
    mpz_urandomm(temp, state,max.get_mpz_t());
    return mpz_class(temp);
}

mpz_class RandomNumberGenerator::rand(int binary_digits){
    mpz_t temp;
    mpz_init(temp);
    mpz_urandomb(temp, state,binary_digits);
    return mpz_class(temp);
}

mpz_class RandomNumberGenerator::generate_prime(unsigned long int n){
    //by the prime number theorem, 10*n attempts should do the trick
    for (int i=0; i<10*n; i++){
        mpz_class temp = rand(n);
        if (mpz_probab_prime_p(temp.get_mpz_t(), MILLER_RABIN_REPEATS) != 0)
            return temp;
    }
    return 1; //failure
}

int legendre_symbol(mpz_class n, mpz_class p){
    //an efficient algorithm using quadratic reciprocity; better than using Euler's criterion
    return jacobi_symbol(n,p);
}

int jacobi_symbol(mpz_class a,mpz_class b){
    if (a % b == 0)
        return 0;
    if (a == 1 || b == 1)
        return 1;
    if (a < 0)
        return jacobi_symbol(-a,b)*((b % 4 == 1)?(1):(-1)); // (-a/b)=(a/b)*(-1)^(b-1 / 4)
    if (a > b)
        return jacobi_symbol(a % b, b); // (a/b) = (a % b / b)

    int count = 0;
    while (a % 2 == 0){
        a /= 2;
        count += 1;
    }
    if (a == 1) // (2^k/b)=(2/b)^k=[(-1)^(b^2-1 / 8)]^k
        if (count % 2 == 0)
            return 1;
        else
            return ((b % 8 == 1 || b % 8 == 7)?(1):(-1));

    int temp = ((a % 4 == 1 || b % 4 == 1)?(1):(-1)); // (a/b)=(b/a)*(-1)^((a-1/4)(b-1/4)) - quadratic reciprocity
    if (count % 2 == 1 && (b % 8 == 3 || b % 8 == 5))
        temp *= -1;
    return temp*jacobi_symbol(b % a, a);
}


//returns x such that x**2 = n. If none exists, returns 0
mpz_class modular_square_root(mpz_class n, mpz_class p){
    if (jacobi_symbol(n,p) != 1)
        return 0;
}
