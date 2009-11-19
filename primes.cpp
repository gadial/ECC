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
    
}
