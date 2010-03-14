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

	print_usg = false;
	do_tests = false;
	pk = false;// string pk_path;
	sk = false;// string sk_path;
	gen_key = false;
	enc = false;// string enc_path;
	dec = false;// string dec_path;
	ec_rand = false;
	ec_prime = false;
	use_ec = false;
	use_predefined_curve = false;
	//string ec_name;


	if (argc <= 1) {
		print_usg = true;
		return;
	}
	for (int i = 1; i < argc; ++i) {
		string str = string(argv[i]);
		if (str == "-enc") {
			if (argc <= ++i) {
				cout << "Specify file to encrypt after -enc !" << endl;
				return;
			}
			enc = true; enc_path = argv[i];
		} else if (str == "-dec") {
			if (argc <= ++i) {
				cout << "Specify file to decrypt after -dec !" << endl;
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
			use_predefined_curve = true;
		} else if (str == "-t") {
			do_tests = true;
			return;
		} else if (str == "-validate") {
			if (argc <= ++i) {
				cout
						<< "Specify an EC declaration after -validate"
						<< endl;
				return;
			}
			validate = true;
			validation_path = argv[i];
		} else if (str == "-use_ec") {
			if (argc <= ++i) {
				cout
						<< "Specify an EC declaration after -validate"
						<< endl;
				return;
			}
			use_ec = true;
			use_ec_path = argv[i];
		}
	}
}

void Cmd::execute() {

	if (print_usg) {
		print_usage();
		return;
	}

	Ellipticcurve* ell = 0;

	// either pk, sk, or keypair-generation should be defined
	if (pk) {
		use_public_key( pk_path);
	} else if (sk) {
		use_private_key( sk_path);
	} else if (gen_key) {
		// now the curve should be defined...
		if (ec_rand) {
			// TODO: implement
		} else if (use_predefined_curve) {
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
				ell = new CurveNISTb163();
			} else if (ec_name == "b283") {
				ell = new CurveNISTb283();
			} else if (ec_name == "b409") {
				ell = new CurveNISTb409();
			} else if (ec_name == "b571") {
				ell = new CurveNISTb571();
			} else {
				cout << "Unknown Curve (" << ec_name << ")!" << endl;
				return;
			}
		} else if (use_ec) {
			ifstream is;
			is.open(use_ec_path.c_str());
			ell = use_elliptic_curve(is);
			is.close();
		}
		elg = new ECC_ElGamal(ell);
		//ell = Ellipticcurve()
		gen_random_keypair();

	} else if (validate) {
		ifstream is;
		is.open(validation_path.c_str());
		ell = use_elliptic_curve(is);
		is.close();
		elg = new ECC_ElGamal(ell);
		if (elg->validate_curve()) {
			cout << "Validation succeeded!" << endl;
		} else {
			cout << "Validation failed!" << endl;
		}
	} else {
		cout << "Either define sk, pk or key-generation" << endl;
	}

	if (enc) {
		encrypt_message( enc_path);
	} else if (dec) {
		decrypt_message( dec_path);
	}

}

void Cmd::print_usage() {
	cout << "Usage:" << endl;
	cout << "ECC {-pk </path/to/pk> -enc </path/to/plaintext> |" << endl <<
			"     -sk </path/to/sk> -dec </path/to/ciphertext> |" << endl <<
			"     -gen_key {-ec_rand | -ec_name 'name' | -use_ec </path/to/ec>}" << endl <<
			"     -validate </path/to/ec>" << endl;
	cout << "Predefined Curves (NIST): [p192|p224|p256|p384|p521|b163|b283|b409|b571]" << endl;
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
		string to_enc = string(str).substr(0, bytes_read);
		//ECC_ElGamal_Ciphertext cipher = elg->encrypt_element(to_enc);
		//os << cipher.to_string(elg->ell) << endl;
		os << elg->encrypt(to_enc) << endl;
		/*
		if (bytes_read < mpl) {
			cipher = elg->encrypt_element();
			os << cipher.to_string(elg->ell) << endl;
		} else {
			cipher = elg->encrypt_element(string(str));
			os << cipher.to_string(elg->ell) << endl;
		}
		*/
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
		string remp = elg->decrypt(line);
		os << remp;
	}
	is.close();
	os.flush();
	os.close();
}

void Cmd::use_public_key(string file_name) {
	ifstream is;
	is.open(file_name.c_str());
	Ellipticcurve* e = use_elliptic_curve(is);
	elg = new ECC_ElGamal(e);
	string line;
	getline(is, line);
	elg->set_public_key(elg->ell->getPointCompressedForm(line));
	is.close();
}

void Cmd::use_private_key(string file_name) {
	ifstream is;
	is.open(file_name.c_str());
	Ellipticcurve* e = use_elliptic_curve(is);
	elg = new ECC_ElGamal(e);
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

	strftime(buffer, 20, "%Y%m%d%H%M", timeinfo);
	string time_str = string(buffer);
	//cout << time_str << endl;

	ofstream of_public;
	of_public.open(("public_key_" + time_str + ".txt").c_str());

	write_elliptic_curve(of_public);

	Coordinate pk = elg->get_public_key();
	of_public << elg->ell->toCompressedForm(pk) << endl;
	of_public.flush();
	of_public.close();

	ofstream of_private;
	of_private.open(("private_key_" + time_str + ".txt").c_str());

	write_elliptic_curve(of_private);

	of_private << elg->get_private_key() << endl;
	of_private.flush();
	of_private.close();
}

Ellipticcurve* Cmd::use_elliptic_curve(ifstream& in) {
	string line;
	// first line: prime/binary
	getline(in, line);
	// TODO: handle binary
	string type = line;

	// 2. line: modulus
	getline(in, line);
	mpz_class modulus;
	modulus.set_str(line, 10);

	// 3. line: a
	getline(in, line);
	mpz_class a;
	a.set_str(line, 10);

	// 4. line: b
	getline(in, line);
	mpz_class b;
	b.set_str(line, 10);

	// 5. line: point on EC
	getline(in, line);
	string compr_point = line;

	// 6. line: EC order
	getline(in, line);
	mpz_class order;
	order.set_str(line, 10);

	Ellipticcurve* ec = 0;
	if (type == "prime") {
		ec = new ECPrime(modulus, a, b);
	} else if (type == "binary") {
		ec = new ECBinary(modulus, a, b);
	} else {
		cout << "EC type neither 'prime' nor 'binary'!" << endl;
	}
	ec->setOrder(order);
	ec->set_point_compressed(compr_point);

	return ec;
	//elg = new ECC_ElGamal(ec);
}

void Cmd::write_elliptic_curve(ofstream& out) {
	Ellipticcurve* ec = elg->ell;
	if (dynamic_cast<ECPrime*>(ec)) {
		out << "prime" << endl;
	} else if (dynamic_cast<ECBinary*>(ec)) {
		out << "binary" << endl;
	} else {
		cerr << "EC has no valid type!" << endl;
	}
	out << ec->mod << endl;
	out << ec->ECC_a << endl;
	out << ec->ECC_b << endl;
	out << ec->toCompressedForm(ec->point) << endl;
	out << ec->getOrder() << endl;
}

Cmd::~Cmd() {
}
