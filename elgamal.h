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

class ECC_ElGamal_Plaintext {
public:
	Coordinate P;
};

class ECC_ElGamal{
public:
    ECC_ElGamal(Ellipticcurve* E);

    static void StringSplit(string str, string delim, vector<string> results) {
		int cutAt;
		while( (cutAt = str.find_first_of(delim)) != str.npos ) {
			if(cutAt > 0) {
				results.push_back(str.substr(0,cutAt));
			}
			str = str.substr(cutAt+1);
		}
		if(str.length() > 0) {
			results.push_back(str);
		}
	}

    Coordinate get_public_key() const{return Q;}
    mpz_class get_private_key() const{return d;}
    void set_keys(Coordinate _Q, mpz_class _d){Q = _Q; d = _d;}
    ECC_ElGamal_Ciphertext encrypt_element(ECC_ElGamal_Plaintext m);
    ECC_ElGamal_Plaintext decrypt_element (ECC_ElGamal_Ciphertext ciphertext);

    string encrypt(string m);
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
     */
    Coordinate to_point(string str);

    /*
     * Returns a point on the EC, with an appropriate padding.
     */
    Coordinate get_point_with_padding(mpz_class str, int padding_length);

private:
    Coordinate Q;
    mpz_class d;
    Ellipticcurve* ell;
    RandomNumberGenerator rand;

    /*
     * Max. length of a point in bytes
     */
    int max_point_length;
};


#endif	/* _ELGAMAL_H */

