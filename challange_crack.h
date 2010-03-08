/* 
 * File:   challange_crack.h
 * Author: gadial
 *
 * Created on February 24, 2010, 11:06 AM
 */

#ifndef _CHALLANGE_CRACK_H
#define	_CHALLANGE_CRACK_H

#include "ecprime.h"

ECPrime candidate_elliptic_curve(int mersenne_prime_number, int D, string ECPoint);

void try_challange_1(ECPrime curve);
void try_challange_2(ECPrime curve);
ECPrime find_suitable_curve();
//void try_all_curves_on_cipher(string ECPoint, string d_string, string cipher);
void try_encryption_and_decryption(ECPrime curve);


void check_example();

#endif	/* _CHALLANGE_CRACK_H */

