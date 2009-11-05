/*
 * coordinates.h
 *
 *  Created on: Nov 5, 2009
 *      Author: bhess
 */

#ifndef COORDINATES_H_
#define COORDINATES_H_

#include <gmpxx.h>
//#include "ellipticcurve.h"

//#define isInfJac(mpz_class m) = m.z

class Jacobian {
public:

	Jacobian() {}
	Jacobian(mpz_class _x, mpz_class _y, mpz_class _z):
		X(_x), Y(_y), Z(_z) {}

	mpz_class X, Y, Z;

};

class Coordinate {
public:

	Coordinate() {}
	Coordinate(mpz_class _x, mpz_class _y):
		X(_x), Y(_y) {}

	mpz_class X, Y;
};

#endif /* COORDINATES_H_ */
