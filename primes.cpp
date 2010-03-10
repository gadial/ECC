#include <gmp.h>
#include <stdlib.h>
#include <time.h>
#include "primes.h"
#include "small_primes.h"
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

zp_int RandomNumberGenerator::generate_modulu_p(mpz_class p){
    return zp_int(rand(p),p);
}

zp_int RandomNumberGenerator::generate_qnr_modulu_p(mpz_class p){
    #define NUMBER_OF_ATTEMPTS_FOR_QNR_GENERATION 20
    zp_int temp;
    for (int i=0; i<NUMBER_OF_ATTEMPTS_FOR_QNR_GENERATION; i++){
        temp = generate_modulu_p(p);
        if (legendre_symbol(temp,p) == -1)
            return temp;
    }
    return 0;
}

mpz_class RandomNumberGenerator::generate_prime_for_discriminant(unsigned long int n, mpz_class D, mpz_class& t, mpz_class& s){
    //find a p such that 4p=t^2+Ds^2 for random t,s
    for (int i=0; i<10*n; i++){
        t = rand_binary_digits(n / 2);
        s = rand_binary_digits(n / 2);
        mpz_class temp = t*t+D*s*s;
        if (temp % 4 == 0){
            mpz_class p = temp / 4;
            if (mpz_probab_prime_p(p.get_mpz_t(), MILLER_RABIN_REPEATS) != 0)
                return p;
        }
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
//    cout << "about to compute jacobi symbol for n = " << n << ", p = " << p << endl;
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
        RandomNumberGenerator gen;
        mpz_class qnr = gen.generate_qnr_modulu_p(p);
        if (qnr == 0) // could not find a QNR
            return -1;
//        cout << "n, p = " << n << ", " << p <<endl;
//        cout << "qnr = " << qnr << endl;
        //now writing p-1 as p-1=2^k*t where t is odd
        int k = 0; //we can safely assume an integer is enough to represent the exponent...
        mpz_class t = p-1;
        while (t % 2 == 0){
            k += 1;
            t /= 2;
        }

        mpz_class z;
        mpz_powm(z.get_mpz_t(),qnr.get_mpz_t(),t.get_mpz_t(),p.get_mpz_t());

//        cout << "k, t = " << k <<", " << t << endl;
//        cout << "z = " << z <<endl;
        //finished the "pre-processing" (up to step 1 in the algorithm, pg. 33)
        mpz_class x,y,b;
        mpz_class temp;
        mpz_powm(y.get_mpz_t(),n.get_mpz_t(),t.get_mpz_t(),p.get_mpz_t());
        temp = (t + 1) / 2;
        mpz_powm(x.get_mpz_t(),n.get_mpz_t(),temp.get_mpz_t(),p.get_mpz_t());
//        cout << "x, y = " << x <<", "<< y << endl;
        mpz_class exp = 1;
        for (int i=0; i<k-2; i++)
            exp *= 2;

        for (int i=0; i<k; i++){
//            cout << "exp = " << exp << endl;
            mpz_powm(b.get_mpz_t(),y.get_mpz_t(),exp.get_mpz_t(),p.get_mpz_t());
//            cout << "i, b = " << i << ", " << b << endl;
            if (b == p-1){
                x = (x*z) % p;
                y = (y*z*z) % p;
            }
            z = (z*z) % p;
            exp /= 2;
//            cout << "x, y, z = " << x << ", " << y <<", " << z << endl;
        }
        return x % p;
}

zp_int modular_square_root(zp_int n){
    //for now, simply use the existing implementation for mpz_class
    return zp_int(modular_square_root(n, n.get_p()),n.get_p());
}
bool is_near_prime(mpz_class p, int smoothness_allowed, mpz_class min_size_allowed){
    int max_prime_num = NUM_SMALL_PRIMES;
    if (smoothness_allowed < NUM_SMALL_PRIMES && smoothness_allowed > 0)
        max_prime_num = smoothness_allowed;
    for (int i=0; i<max_prime_num; i++){
        if (p % small_primes[i] == 0)
            p /= small_primes[i];
    }

    if (p<min_size_allowed)
        return false;
    return mpz_probab_prime_p(p.get_mpz_t(), MILLER_RABIN_REPEATS);
}

mpz_class mersenne_prime(int n){
    #define MERSENNE_NUMERS_SIZE 20
    int p_values[MERSENNE_NUMERS_SIZE] = {2,3,5,7,13,17,19,31,61,89,107,127,521,607,1279,2203,2281,3217,4253,4423};
    if (n > MERSENNE_NUMERS_SIZE)
        return 0;
    mpz_class result;
    mpz_class base = 2;
    mpz_pow_ui(result.get_mpz_t(),base.get_mpz_t(),p_values[n-1]);
    return (result - 1);
}

//finds t,s such that t^2+|D|y^2=4p, D is negative and equivalent to 0 or 1 modulo 4, and |D| < 4p
bool extended_cornacchia(mpz_class p, int D, mpz_class& t,mpz_class& s){
    if (p == 2){
        switch (D+8){
            case 0: t = 0; s = 1; return true;
            case 1: t = 1; s = 1; return true;
            case 4: t = 2; s = 1; return true;
            default: return false;
        }
    }

    if (legendre_symbol(D,p) == -1)
        return false;

    mpz_class d = (D < 0)?(-D):(D);
    mpz_class x0 = modular_square_root(D,p);
    if (x0 % 2 != d % 2)
        x0 = p - x0;

    mpz_class a,b,l,r;
    a = 2*p;
    b = x0;
    mpz_root(l.get_mpz_t(),p.get_mpz_t(),2);
    l*=2;

    while (b > l){
        r = a % b;
        a = b;
        b = r;
    }

    
    mpz_class c = 4*p-b*b;

    if (c % d != 0)
        return false;
    c/= d;

    mpz_class result;
    if (mpz_root(result.get_mpz_t(), c.get_mpz_t(),2)){
        t = b;
        s = result;
        return true;
    }
    return false;
}

bool is_near_prime_by_min_max(mpz_class p, mpz_class first_divisor_min, int second_divisor_max){
    for (int i=2; i<second_divisor_max; i++){
        while (p % i == 0)
            p /= i;
    }
    return ((p > first_divisor_min) && mpz_probab_prime_p(p.get_mpz_t(), MILLER_RABIN_REPEATS));
}
