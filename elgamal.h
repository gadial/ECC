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
    ECC_ElGamal(Ellipticcurve E);

    Coordinate get_public_key() const{return Q;}
    mpz_class get_private_key() const{return d;}
    void set_keys(Coordinate _Q, mpz_class _d){Q = _Q; d = _d;}
    ECC_ElGamal_Ciphertext ecnrypt_element(string m) const;
    string decrypt_element (ECC_ElGamal_Ciphertext ciphertext) const;

    string encrypt(string m) const;
    string decrypt(string c) const;

private:
    Coordinate Q;
    mpz_class d;
};


#endif	/* _ELGAMAL_H */

