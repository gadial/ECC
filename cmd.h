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

	void encrypt_message(string file_name);
	void decrypt_message(string file_name);
	void gen_random_keypair();
	void use_public_key(string file_name);
	void use_private_key(string file_name);

	void print_usage();

	virtual ~Cmd();

	ECC_ElGamal* elg;
};

#endif /* CMD_H_ */
