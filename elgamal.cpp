/*
 * elgamal.cpp
 *
 *  Created on: Jan 27, 2010
 *      Author: bhess
 */

#include "elgamal.h"

ECC_ElGamal::ECC_ElGamal(Ellipticcurve* E) {}

ECC_ElGamal_Ciphertext ECC_ElGamal::ecnrypt_element(ECC_ElGamal_Plaintext m) {
	ECC_ElGamal_Ciphertext eco;
	Coordinate M = ell->getPoint(m.msg_mpz, false);
	Coordinate P = ell->point;
	mpz_class n = ell->getOrder();
	mpz_class k = rand.rand(n - 1);
	eco.C1 = ell->pointMultiplication(P, k);
	Coordinate kQ = ell->pointMultiplication(Q, k);
	eco.C2 = ell->addition(M, kQ);
	return eco;
}

ECC_ElGamal_Plaintext ECC_ElGamal::decrypt_element (ECC_ElGamal_Ciphertext ciphertext) {
	ECC_ElGamal_Plaintext ep;
	Coordinate dC1 = ell->pointMultiplication(ciphertext.C1, d);
	ep.msg_mpz = ell->subtraction(ciphertext.C2, dC1).X;
	return ep;
}
