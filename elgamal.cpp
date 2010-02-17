/*
 * elgamal.cpp
 *
 *  Created on: Jan 27, 2010
 *      Author: bhess
 */

#include "elgamal.h"

ECC_ElGamal::ECC_ElGamal(Ellipticcurve* E) { ell = E;}

ECC_ElGamal_Ciphertext ECC_ElGamal::encrypt_element(ECC_ElGamal_Plaintext m) {
	//cout << "Message: " << m.P << endl;
	ECC_ElGamal_Ciphertext eco;
	Coordinate M = m.P;
	Coordinate P = ell->point;
	mpz_class n = ell->getOrder();
	mpz_class k = rand.rand(n - 1);
	eco.C1 = ell->pointMultiplication(P, k);
	Coordinate kQ = ell->pointMultiplication(Q, k);
	eco.C2 = ell->addition(M, kQ);
	//Coordinate ccc = ell->getPoint(M.X);
	//cout << "getpoint c2: " << M << endl;
	return eco;
}

ECC_ElGamal_Ciphertext ECC_ElGamal::encrypt_element(string m) {
	return encrypt_element(ECC_ElGamal_Plaintext(to_point(m)));
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

		//cout << "C1: " << cipher.C1 << endl << "C1-comp: " << cipher.C1.toCompressedForm() << endl;
		//cout << "C2: " << cipher.C2 << endl << "C2-comp: " << cipher.C2.toCompressedForm() << endl;
		/*
		res.append(cipher.C1.toCompressedForm());
		res.append(",");
		res.append(cipher.C2.toCompressedForm());
		res.append(",");
		*/
		res.append(cipher.to_string());
	}
	return res;
}

string ECC_ElGamal::decrypt(string c) {

	vector<string> ciphers = StringSplit(c, ",");
	//cout << "size: " << ciphers.size() << endl;
	string res;
	for (unsigned int i = 0; i < ciphers.size(); i += 2) {
		ECC_ElGamal_Ciphertext ec;
		ec.C1 = ell->getPointCompressedForm(ciphers[i]);
		ec.C2 = ell->getPointCompressedForm(ciphers[i + 1]);
		//cout << "C1: " << ciphers[i] << endl << ec.C1 << endl;
		//cout << "C2: " << ciphers[i + 1] << endl << ec.C2 << endl;
		ECC_ElGamal_Plaintext ep = decrypt_element(ec);
		//cout << ep.P << endl;
		// removing padding...

		/*
		ep.P.X >>= 8;
		res.append(to_string(ep.P.X));
		*/

		res.append(remove_padding(ep));
	}
	return res;
}

string ECC_ElGamal::remove_padding(ECC_ElGamal_Plaintext& ep) {
	return to_string(ep.P.X >> 8);
}

vector<ECC_ElGamal_Plaintext> ECC_ElGamal::split_msg(string msg) {
	int str_length_bits = msg.length();
	int max_point_length_effective = get_max_point_length();
	//cout << max_point_length_effective << endl;
	int vec_el = (str_length_bits % max_point_length_effective == 0 ?
			str_length_bits / max_point_length_effective : (str_length_bits / max_point_length_effective) + 1);
	vector<ECC_ElGamal_Plaintext> res(vec_el);

	for (int i = 0; i < vec_el; ++i) {
		string s = msg.substr(i * max_point_length_effective, max_point_length_effective);
		ECC_ElGamal_Plaintext ep;
		//cout << "Part " << i << ":" << s << endl;
		ep.P = to_point(s);
		res[i] = ep;
	}
	return res;
}

int ECC_ElGamal::get_max_point_length() {
	return (ell->get_bits() / 8) - 2;
}

Coordinate ECC_ElGamal::to_point(string str) {
	mpz_class res;
	// 8 bit encoding of a character, as in ASCII...
	//res.set_str(str, 256);
	for (int i = str.length() - 1; i >= 0; --i) {
		res |= (int)str[i];
		res <<= 8;
	}
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
	//cout << "Str: " << to_string(str_copy) << endl;
	return c;
}

void ECC_ElGamal::generate_random_keypair() {
	d = rand.rand(ell->getOrder() - 1);
	Q = ell->pointMultiplication(ell->point, d);
}

string ECC_ElGamal::to_string(mpz_class mpz) {
	string res;
	while (mpz != 0) {
		res += (char)(mpz.get_ui() & 0xFF);
		mpz >>= 8;
	}
	return res;
}
