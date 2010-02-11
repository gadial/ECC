/*
 * elgamal.cpp
 *
 *  Created on: Jan 27, 2010
 *      Author: bhess
 */

#include "elgamal.h"

ECC_ElGamal::ECC_ElGamal(Ellipticcurve* E) { ell = E;}

ECC_ElGamal_Ciphertext ECC_ElGamal::encrypt_element(ECC_ElGamal_Plaintext m) {
	ECC_ElGamal_Ciphertext eco;
	Coordinate M = m.P;
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
	ep.P = ell->subtraction(ciphertext.C2, dC1);
	return ep;
}

string ECC_ElGamal::encrypt(string m) {
	vector<ECC_ElGamal_Plaintext> co = split_msg(m);
	string res;
	for (unsigned int i = 0; i < co.size(); ++i) {
		ECC_ElGamal_Ciphertext cipher = encrypt_element(co[i]);
		/*
		 * TODO: some better output format?
		 */
		res.append(cipher.C1.toCompressedForm());
		res.append(",");
		res.append(cipher.C2.toCompressedForm());
		res.append(",");
	}
	return res;
}

string ECC_ElGamal::decrypt(string c) {
	vector<string> ciphers;
	StringSplit(c, ",", ciphers);
	string res;
	for (unsigned int i = 0; i < ciphers.size(); i += 2) {
		ECC_ElGamal_Ciphertext ec;
		ec.C1 = Coordinate::fromCompressedForm(ciphers[i]);
		ec.C2 = Coordinate::fromCompressedForm(ciphers[i + 1]);
		ECC_ElGamal_Plaintext ep = decrypt_element(ec);
		// removing padding...
		ep.P.X >>= 8;
		res.append(ep.P.X.get_str(256));
	}
	return res;
}

vector<ECC_ElGamal_Plaintext> ECC_ElGamal::split_msg(string msg) {
	int str_length_bits = msg.length();
	int max_point_length_effective = str_length_bits - 1;
	int vec_el = (str_length_bits % max_point_length_effective == 0 ?
			str_length_bits / max_point_length_effective : (str_length_bits / max_point_length_effective) + 1);
	vector<ECC_ElGamal_Plaintext> res(vec_el);

	for (int i = 0; i < vec_el; ++i) {
		string s = msg.substr(i * max_point_length_effective, max_point_length_effective);
		ECC_ElGamal_Plaintext ep;
		ep.P = to_point(s);
	}
	return res;
}

Coordinate ECC_ElGamal::to_point(string str) {
	mpz_class res;
	// 8 bit encoding of a character, as in ASCII...
	res.set_str(str, 256);
	return get_point_with_padding(res, 1);
}

Coordinate ECC_ElGamal::get_point_with_padding(mpz_class str, int padding_length) {
	mpz_class str_copy;
	Coordinate c;
	do {
		// setting random padding until there exists a point
		str_copy = str;
		mpz_class pad = rand.rand(255);
		str_copy <<= 8;
		str_copy |= pad;
		c = ell->getPoint(str_copy, false);
	} while(c.isInfinite());
	return c;
}
