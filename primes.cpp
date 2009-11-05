#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "primes.h"

#include <stdio.h>
#define MILLER_RABIN_REPEATS 5

mpz_class generate_prime(unsigned long int n){
    gmp_randstate_t random_state;
    gmp_randinit_mt(random_state);
    gmp_randseed_ui(random_state,time(NULL)); //for now, use a random seed
    mpz_t temp;
    mpz_init(temp);

    //by the prime number theorem, 10*n attempts should do the trick
    for (int i=0; i<10*n; i++){
        mpz_urandomb(temp, random_state,n);
        if (mpz_probab_prime_p(temp, MILLER_RABIN_REPEATS) != 0)
            return mpz_class(temp);
    }
    return 1; //failure
}