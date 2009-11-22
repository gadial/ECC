/*
 * coordinates.cpp
 *
 *  Created on: Nov 12, 2009
 *      Author: bhess
 */

#include "coordinates.h"


Coordinate::Coordinate(const Jacobian& jac, const mpz_class mod) {
	mpz_class const_2 = 2, const_3 = 3;

	mpz_powm(X.get_mpz_t(), jac.Z.get_mpz_t(), const_2.get_mpz_t(), mod.get_mpz_t());
	mpz_powm(Y.get_mpz_t(), jac.Z.get_mpz_t(), const_3.get_mpz_t(), mod.get_mpz_t());
	mpz_invert(X.get_mpz_t(), X.get_mpz_t(), mod.get_mpz_t());
	mpz_invert(Y.get_mpz_t(), Y.get_mpz_t(), mod.get_mpz_t());
	mpz_mul(X.get_mpz_t(), jac.X.get_mpz_t(), X.get_mpz_t());
	mpz_mul(Y.get_mpz_t(), jac.Y.get_mpz_t(), Y.get_mpz_t());

	mpz_mod(X.get_mpz_t(), X.get_mpz_t(), mod.get_mpz_t());
	mpz_mod(Y.get_mpz_t(), Y.get_mpz_t(), mod.get_mpz_t());
}

Coordinate::Coordinate(const LD& ld, const mpz_class mod) {
	mpz_class const_2 = 2;
	//mpz_powm(X.get_mpz_t(), jac.Z.get_mpz_t(), const_2.get_mpz_t(), mod.get_mpz_t());
	mpz_powm(Y.get_mpz_t(), ld.Z.get_mpz_t(), const_2.get_mpz_t(), mod.get_mpz_t());
	mpz_invert(X.get_mpz_t(), ld.Z.get_mpz_t(), mod.get_mpz_t());
	mpz_invert(Y.get_mpz_t(), Y.get_mpz_t(), mod.get_mpz_t());
	mpz_mul(X.get_mpz_t(), ld.X.get_mpz_t(), X.get_mpz_t());
	mpz_mul(Y.get_mpz_t(), ld.Y.get_mpz_t(), Y.get_mpz_t());

	mpz_mod(X.get_mpz_t(), X.get_mpz_t(), mod.get_mpz_t());
	mpz_mod(Y.get_mpz_t(), Y.get_mpz_t(), mod.get_mpz_t());
}

Coordinate::Coordinate(const char* _x, int basex,
			const char* _y, int basey) {
	X.set_str(_x, basex);
	Y.set_str(_y, basey);
}
