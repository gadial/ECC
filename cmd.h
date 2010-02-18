/*
 * cmd.h
 *
 *  Created on: Feb 17, 2010
 *      Author: bhess
 *
 *      This class represents the command line options
 *      for the Elliptic curve program
 */

#ifndef CMD_H_
#define CMD_H_

using namespace std;

#include "elgamal.h"
#include <iostream>

class Cmd {
public:
	Cmd(int argc, char** argv);

	/*
	 * Encrypts the message given by a file path 'file_name'
	 * Encodes the message to points on the EC
	 *
	 * Output: ciphertext-points (in compressed format), line by line
	 * as 'file_name'.enc
	 */
	void encrypt_message(string file_name);

	/*
	 * Decrypts a given file with ciphertext-points
	 *
	 * Output: the original plaintext as 'file_name'.dec
	 */
	void decrypt_message(string file_name);

	/*
	 * Generates a random keypair.
	 *
	 * The generated files include the EC definition
	 *
	 * Output (public_key_'time'.txt, private_key_'time'.txt)
	 *
	 * 1. line: prime/binary curve
	 * 2. line: Elliptic curve modulus
	 * 3. line: parameter 'a' of elliptic curve
	 * 4. line: parameter 'b' of elliptic curve
	 * 5. line: point on elliptic curve (compressed format)
	 * 6. line: EC order
	 * 7. line: public key / private key
	 */
	void gen_random_keypair();

	/*
	 * Using an existing public key file
	 */
	void use_public_key(string file_name);

	/*
	 * Using an existing private key file
	 */
	void use_private_key(string file_name);

	void print_usage();
	void execute();

	bool do_tests;

	virtual ~Cmd();

	ECC_ElGamal* elg;

private:

	/*
	 * using the definition of an elliptic curve
	 *
	 * initializes elg
	 */
	void use_elliptic_curve(ifstream& in);

	/*
	 * writes the definition of the elliptic curve to the out-stream
	 */
	void write_elliptic_curve(ofstream& out);

	bool pk;
	string pk_path;
	bool sk;
	string sk_path;
	bool gen_key;
	bool enc;
	string enc_path;
	bool dec;
	string dec_path;
	bool ec_rand;
	bool ec_prime;
	string ec_name;
	bool print_usg;
};

#endif /* CMD_H_ */
