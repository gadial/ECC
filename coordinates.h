/*
 * coordinates.h
 *
 *  Created on: Nov 5, 2009
 *      Author: bhess
 */

#ifndef COORDINATES_H_
#define COORDINATES_H_

#include <gmpxx.h>

class Jacobian {
public:

	Jacobian() {}
	Jacobian(mpz_class _x, mpz_class _y, mpz_class _z):
		X(_x), Y(_y), Z(_z) {}

	mpz_class X, Y, Z;
};

#endif /* COORDINATES_H_ */
