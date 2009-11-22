#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "primes.h"

#include <stdio.h>
#include <iostream>

#define MILLER_RABIN_REPEATS 5

using namespace std;

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

mpz_class RandomNumberGenerator::rand_binary_digits(int binary_digits){
    mpz_t temp;
    mpz_init(temp);
    mpz_urandomb(temp, state,binary_digits);
    return mpz_class(temp);
}

mpz_class RandomNumberGenerator::generate_prime(unsigned long int n){
    //by the prime number theorem, 10*n attempts should do the trick
    for (int i=0; i<10*n; i++){
        mpz_class temp = rand_binary_digits(n);
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


//returns x such that x**2 = n. If none exists, returns 0. On failure, returns -1
mpz_class modular_square_root(mpz_class n, mpz_class p){ // we follow Cohen's computational number theory book, pg. 32
    if (jacobi_symbol(n,p) != 1)
        return 0;
    mpz_class result;
    if (p % 4 == 3){ // a very simple case
        mpz_class exp = (p+1) / 4;
        mpz_powm(result.get_mpz_t(),n.get_mpz_t(),exp.get_mpz_t(),p.get_mpz_t());
        return result;
    }
    if (p % 8 == 5){ // still a simple case
        mpz_class exp = (p-1) / 4;
        mpz_powm(result.get_mpz_t(),n.get_mpz_t(),exp.get_mpz_t(),p.get_mpz_t());
        if (result == 1){
//            cout << "case 1" << endl;
            exp = (p+3) / 8;
            mpz_powm(result.get_mpz_t(),n.get_mpz_t(),exp.get_mpz_t(),p.get_mpz_t());
            return result;
        }
        else{
//            cout << "case 2" << endl;
            exp = (p - 5) / 8;
            mpz_class temp_n = n*4;
            mpz_powm(result.get_mpz_t(),temp_n.get_mpz_t(),exp.get_mpz_t(),p.get_mpz_t());
            return ((2*n*result) % p);
        }
    }
        //we now remain in the "hard" case of p % 8 == 1, and use Shanks-Tonelli

        //first step - obtain a non-quadratic residue
        int number_of_attempts = 20;
        mpz_class qnr = 0;
        mpz_class temp;
        RandomNumberGenerator gen;
        for (int i=0; i<number_of_attempts; i++){
            temp = gen.rand(p-2) + 2;
            if (legendre_symbol(temp,p) == -1){
                qnr = temp;
                break;
            }
        }
        if (qnr == 0) // could not find a QNR
            return -1;

        //now writing p-1 as p-1=2^k*t where t is odd
        int k = 0; //we can safely assume an integer is enough to represent the exponent...
        mpz_class t = p-1;
        while (t % 2 == 0){
            k += 1;
            t /= 2;
        }

        mpz_class z;
        mpz_powm(z.get_mpz_t(),qnr.get_mpz_t(),t.get_mpz_t(),p.get_mpz_t());

        //finished the "pre-processing" (up to step 1 in the algorithm, pg. 33)
        mpz_class x,y,b;
        mpz_powm(y.get_mpz_t(),n.get_mpz_t(),t.get_mpz_t(),p.get_mpz_t());
        temp = (t + 1) / 2;
        mpz_powm(x.get_mpz_t(),n.get_mpz_t(),temp.get_mpz_t(),p.get_mpz_t());

        mpz_class exp;
        for (int i=0; i<k-2; i++)
            exp *= 2;

        for (int i=0; i<k; i++){         
            mpz_powm(b.get_mpz_t(),y.get_mpz_t(),exp.get_mpz_t(),p.get_mpz_t());
            if (b == p-1){
                x = (x*z) % p;
                y = (y*z*z) % p;
            }
            z = (z*z) % p;
            exp /= 2;
        }
        return x % p;
}
