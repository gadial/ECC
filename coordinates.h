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
	Coordinate(const char* _x, int basex,
			const char* _y, int basey);
	Coordinate(mpz_class _x, mpz_class _y):
		X(_x), Y(_y) {}
	Coordinate(const Jacobian& jac, const mpz_class mod);

        //returns the point at infinity, as is represented by this class in the context of elliptic curves
        static Coordinate infinity(){return Coordinate(0,0);}

	bool operator==(const Coordinate& eqTo) {
		return X == eqTo.X && Y == eqTo.Y;
	}

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
