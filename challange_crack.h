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

void try_challange_1();
    
void try_all_curves_on_cipher(string ECPoint, string d_string, string cipher);
#endif	/* _CHALLANGE_CRACK_H */

