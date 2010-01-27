/* 
 * File:   elgamal.h
 * Author: gadial
 *
 * Created on January 27, 2010, 2:51 PM
 */

#ifndef _ELGAMAL_H
#define	_ELGAMAL_H
#include "ellipticcurve.h"
class ECC_ElGamal_Ciphertext{
public:
    Coordinate C1, C2;
};

class ECC_ElGamal{
public:
    ECC_ElGamal(mpz_class p, Ellipticcurve E, Coordinate P, mpz_class n);

    Coordinate get_public_key() const{return Q;}
    mpz_class get_private_key() const{return d;}
    void set_keys(Coordinate _Q, mpz_class _d){Q = _Q; d = _d;}
    ECC_ElGamal_Ciphertext ecnrypt(string m) const;
    string decrypt (ECC_ElGamal_Ciphertext ciphertext, mpz_class d) const;

private:
    Coordinate Q;
    mpz_class d;
};


#endif	/* _ELGAMAL_H */

