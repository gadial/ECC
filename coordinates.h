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

class Coordinate;
class Jacobian;

class Coordinate {
public:

	Coordinate() {}
	Coordinate(mpz_class _x, mpz_class _y):
		X(_x), Y(_y) {}
	Coordinate(const Jacobian& jac, const mpz_class mod);

    bool isInfinite() {
    	return X == 0 && Y == 0;
    }

	mpz_class X, Y;
};

class Jacobian {
public:

	Jacobian() {}
	Jacobian(mpz_class _x, mpz_class _y, mpz_class _z):
		X(_x), Y(_y), Z(_z) {}
        Jacobian(const Coordinate& rhs): //TODO: is this really the correct conversion?
                X(rhs.X), Y(rhs.Y), Z(1) {}


    bool isInfinite() {
    	return X == 1 && Y == 1 && Z == 0;
    }

	mpz_class X, Y, Z;

};

#endif /* COORDINATES_H_ */
