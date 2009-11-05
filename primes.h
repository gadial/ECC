/* 
 * File:   primes.h
 * Author: gadial
 *
 * Created on November 5, 2009, 2:47 PM
 */

#ifndef _PRIMES_H
#define	_PRIMES_H
#include <gmpxx.h>

//generates a random prime not larger than 2^n-1
mpz_class generate_prime(unsigned long int n);


#endif	/* _PRIMES_H */

