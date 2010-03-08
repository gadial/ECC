/* 
 * File:   elgamal.h
 * Author: gadial
 *
 * Created on January 27, 2010, 2:51 PM
 */

#ifndef _ELGAMAL_H
#define	_ELGAMAL_H
#include "ellipticcurve.h"

static vector<string> StringSplit(string str, string delim) {
	int cutAt;
	vector<string> results;
	while ((cutAt = str.find_first_of(delim)) != str.npos) {
		if (cutAt > 0) {
			results.push_back(str.substr(0, cutAt));
		}
		str = str.substr(cutAt + 1);
	}
	if (str.length() > 0) {
		results.push_back(str);
	}
	return results;
}

class ECC_ElGamal_Ciphertext{
public:
	static ECC_ElGamal_Ciphertext from_string(string str, Ellipticcurve* el) {
		ECC_ElGamal_Ciphertext res;
		//cout << str << endl;
		vector<string> v = StringSplit(str, ",");
		res.C1 = el->getPointCompressedForm(v[0]);
		res.C2 = el->getPointCompressedForm(v[1]);
		return res;
	}

    Coordinate C1, C2;

    string to_string() {
    	string res;
    	res.append(C1.toCompressedForm());
    	res.append(",");
    	res.append(C2.toCompressedForm());
    	//res.append();
    	return res;
    }
};

class ECC_ElGamal_Plaintext {
public:
	ECC_ElGamal_Plaintext() {}
	ECC_ElGamal_Plaintext(Coordinate _P) : P(_P) {}

	Coordinate P;
};

class ECC_ElGamal{
public:
    ECC_ElGamal(Ellipticcurve* E);

    Coordinate get_public_key() const{return Q;}
    mpz_class get_private_key() const{return d;}
    void set_keys(Coordinate _Q, mpz_class _d){Q = _Q; d = _d;}
    void set_public_key(Coordinate _Q) {Q = _Q;}
    void set_private_key(mpz_class _d) {d = _d;}
    void set_keys_from_private_key(mpz_class _d){d = _d; Q = ell->pointMultiplication(ell->point, d);}

    /*
     * Generates a random pk/sk key pair according to Algorithm 1.12 in
     * "Guide to Elliptic Curve Cryptography"
     */
    void generate_random_keypair();
    ECC_ElGamal_Ciphertext encrypt_element(ECC_ElGamal_Plaintext m);
    ECC_ElGamal_Ciphertext encrypt_element(string m);

    ECC_ElGamal_Plaintext decrypt_element (ECC_ElGamal_Ciphertext ciphertext);

    /*
     * Encrypts a given string 'm'
     * 'm' must be no longer than 'get_max_point_length()'
     */
    string encrypt(string m);

    /*
     * decrypts 'c'
     * 'c' is a ciphertext in the following format: C1,C2
     * Where C1 and C2 are in compressed format
     */
    string decrypt(string c);

    /*
     * Splits the message string into coordinates on the EC
     * integers, with 1 byte random padding (p)
     *     str_1   p      str_2   p
     * [|---------|-|],[---------|-|],...
     *
     */
    vector<ECC_ElGamal_Plaintext> split_msg(string msg);

    /*
     * Converts a string of length 'max_point_length'-1 to an gmp integer
     *
     * More specifically:
     *
     * Encodes a string to a point on the Elliptic Curve, including
	 * one byte random padding. In the following way:
	 *
	 * The string corresponds to the X-coordinate of the point,
	 * random padding is choosen in such a way that a Y coordinate exists
	 *
	 * --------------------------
	 * Format of the X-Coordinate
	 * --------------------------
	 *
	 * ASCII(char) means the ascii code of a character (in binary format)
	 * rand() is a random (one byte) character
	 *
	 * Example: "abc" <-> ASCII(c)|ASCII(b)|ASCII(a)|ASCII(rand())
	 */
    Coordinate to_point(string str);

    /*
     * Returns a point on the EC, with an appropriate padding.
     */
    Coordinate get_point_with_padding(mpz_class str, int padding_length);

    void print_keypair() {
    	cout << "private key: " << get_private_key() << endl;
    	cout << "public key: " << get_public_key().X << "," << get_public_key().Y << endl;
    }

    string remove_padding(ECC_ElGamal_Plaintext& ep);

    int get_max_point_length();
    Ellipticcurve* ell;

private:
    Coordinate Q;
    mpz_class d;
    RandomNumberGenerator rand;

    string to_string(mpz_class mpz);

    /*
     * Max. length of a point in bytes
     */
    int max_point_length;
};


#endif	/* _ELGAMAL_H */

