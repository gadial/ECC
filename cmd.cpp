/*
 * cmd.cpp
 *
 *  Created on: Feb 17, 2010
 *      Author: bhess
 */

#include "cmd.h"
#include "curvesnist.h"

#include <fstream>

Cmd::Cmd(int argc, char** argv) {

	bool pk = false; string pk_path;
	bool sk = false; string sk_path;
	bool gen_key = false;
	bool enc = false; string enc_path;
	bool dec = false; string dec_path;
	bool ec_rand = false; bool ec_prime = false;
	string ec_name;


	if (argc <= 1) {
		print_usage();
		return;
	}
	for (int i = 1; i < argc; ++i) {
		string str = string(argv[i]);
		if (str == "-e") {
			if (argc <= ++i) {
				cout << "Specify file to encrypt after -e !" << endl;
				return;
			}
			enc = true; enc_path = argv[i];
		} else if (str == "-d") {
			if (argc <= ++i) {
				cout << "Specify file to decrypt after -d !" << endl;
				return;
			}
			dec = true; dec_path = argv[i];
		} else if (str == "-pk") {
			if (argc <= ++i) {
				cout << "Specify public key after -pk !" << endl;
				return;
			}
			pk = true;
			pk_path = argv[i];
		} else if (str == "-sk") {
			if (argc <= ++i) {
				cout << "Specify private key after -sk !" << endl;
				return;
			}
			sk = true;
			sk_path = argv[i];
		} else if (str == "-gen_key") {
			gen_key = true;
		} else if (str == "-ec_rand") {
			if (argc <= ++i) {
				cout << "Specify type of random curve after -ec_rand: either prime or binary" << endl;
				return;
			}
			ec_rand = true;
			ec_prime = (string(argv[i]) == "prime");
		} else if (str == "-ec_name") {
			if (argc <= ++i) {
				cout
						<< "Specify which random curve after -ec_name"
						<< endl;
				return;
			}
			ec_name = argv[i];
		}
	}

	Ellipticcurve* ell;

	// either pk, sk, or keypair-generation should be defined
	if (pk) {
		use_public_key(pk_path);
		// TODO: load elliptic curve from config...
	} else if (sk) {
		use_private_key(sk_path);
		// TODO: load...
	} else if (gen_key) {
		// now the curve should be defined...
		if (ec_rand) {
			// TODO: implement
		} else {
			if (ec_name == "p192") {
				ell = new CurveNISTp192();
			} else if (ec_name == "p224") {
				ell = new CurveNISTp224();
			} else if (ec_name == "p256") {
				ell = new CurveNISTp256();
			} else if (ec_name == "p384") {
				ell = new CurveNISTp384();
			} else if (ec_name == "p521") {
				ell = new CurveNISTp521();
			} else if (ec_name == "b163") {
				ell = new CurveNISTb163;
			} else {
				cout << "Unknown Curve (" << ec_name << ")!" << endl;
				return;
			}
		}
		elg = new ECC_ElGamal(ell);
		//ell = Ellipticcurve()
		gen_random_keypair();

	} else {
		cout << "Either define sk, pk or key-generation" << endl;
	}


	if (enc) {
		encrypt_message(enc_path);
	} else if (dec) {
		decrypt_message(dec_path);
	}

}

void Cmd::print_usage() {
	cout << "Usage:" << endl;
	cout << "ECC {-pk </path/to/pk> -enc </path/to/plaintext> |" << endl <<
			"     -sk </path/to/sk> -dec </path/to/ciphertext> |" << endl <<
			"     -gen_key {-ec_rand | -ec_name 'name'} [-enc ... | -dec ...]}" << endl << endl;
	cout << "Predefined Curves (NIST): [p192|p224|p256|p384|p521|b163]" << endl;
}

void Cmd::encrypt_message(string file_name) {
	//ECC_ElGamal* elg = new ECC_ElGamal(new CurveNISTp192());
	int mpl = elg->get_max_point_length();
	char str[mpl];

	ifstream is;
	is.open(file_name.c_str(), ios::binary);
	ofstream os;
	os.open((file_name + ".enc").c_str());

	//is.seekg (0, ios::end);
	//int length = is.tellg();
	//is.seekg (0, ios::beg);

	while (is.good()) {
		is.read(str, mpl);
		streamsize bytes_read = is.gcount();
		ECC_ElGamal_Ciphertext cipher = elg->encrypt_element(
				bytes_read < mpl ? string(str).substr(0, bytes_read) : string(
						str));
		os << cipher.to_string() << endl;
	}
	is.close();
	os.flush();
	os.close();
}

void Cmd::decrypt_message(string file_name) {
	ifstream is;
	is.open(file_name.c_str(), ios::binary);
	ofstream os;
	os.open((file_name + ".dec").c_str());
	string line;
	while (!is.eof()) {
		getline(is, line);
		if (line.empty()) continue;
		ECC_ElGamal_Plaintext ep = elg->decrypt_element(ECC_ElGamal_Ciphertext::from_string(line, elg->ell));
		string remp = elg->remove_padding(ep);
		os << remp;
	}
	is.close();
	os.flush();
	os.close();
}

void Cmd::use_public_key(string file_name) {
	ifstream is;
	is.open(file_name.c_str());
	string line;
	getline(is, line);
	elg->set_public_key(elg->ell->getPointCompressedForm(line));
	is.close();
}

void Cmd::use_private_key(string file_name) {
	ifstream is;
	is.open(file_name.c_str());
	string line;
	getline(is, line);
	mpz_class sk;
	sk.set_str(line, 10);
	elg->set_private_key(sk);
	is.close();
}

void Cmd::gen_random_keypair() {
	elg->generate_random_keypair();

	time_t rawtime;
	struct tm * timeinfo;
	char buffer[20];

	time(&rawtime);
	timeinfo = localtime(&rawtime);

	strftime(buffer, 20, "%Y%m%d%H%M%S", timeinfo);
	string time_str = string(buffer);
	//cout << time_str << endl;

	ofstream of_public;
	of_public.open(("public_key_" + time_str + ".txt").c_str());
	Coordinate pk = elg->get_public_key();
	of_public << pk.toCompressedForm() << endl;
	of_public.flush();
	of_public.close();

	ofstream of_private;
	of_private.open(("private_key_" + time_str + ".txt").c_str());
	of_private << elg->get_private_key() << endl;
	of_private.flush();
	of_private.close();
}

Cmd::~Cmd() {
}
