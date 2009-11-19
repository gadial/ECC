/* 
 * File:   primes.h
 * Author: gadial
 *
 * Created on November 5, 2009, 2:47 PM
 */

#ifndef _PRIMES_H
#define	_PRIMES_H
#include <gmpxx.h>
class RandomNumberGenerator{ //naive. very naive
public:
    RandomNumberGenerator(unsigned long seed = 0);
    mpz_class rand(mpz_class max);
    mpz_class rand_binary_digits(int binary_digits);
    //generates a random prime not larger than 2^n-1
    mpz_class generate_prime(unsigned long int n);
private:
    gmp_randstate_t state;
};


//returns 1 if n is a quadratic residue modulo a prime p, -1 if it is not, and 0 if p divides n
int legendre_symbol(mpz_class n, mpz_class p);

//returns 1 if a is a quadratic residue modulo an odd b, -1 if it is not, and 0 if b divides a
int jacobi_symbol(mpz_class a,mpz_class b);
   
//returns x such that x**2 = n. If none exists, returns 0
mpz_class modular_square_root(mpz_class n, mpz_class p);

#endif	/* _PRIMES_H */

