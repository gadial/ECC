/*
 * coordinates.h
 *
 *  Created on: Nov 5, 2009
 *      Author: bhess
 */

#ifndef COORDINATES_H_
#define COORDINATES_H_
#include <iostream>

#include <gmpxx.h>
using namespace std;
//#include "ellipticcurve.h"
//#define isInfJac(mpz_class m) = m.z

class Coordinate;
class Jacobian;
class LD;

class Coordinate {
public:

	static Coordinate fromCompressedForm(string comp) {
		// Todo: implement
	}

	Coordinate() {}
	Coordinate(const char* _x, int basex,
			const char* _y, int basey);
	Coordinate(mpz_class _x, mpz_class _y):
		X(_x), Y(_y) {}
	Coordinate(const Jacobian& jac, const mpz_class mod);
	Coordinate(const LD& ld, const mpz_class mod);

        //returns the point at infinity, as is represented by this class in the context of elliptic curves
        static Coordinate infinity(){return Coordinate(0,0);}

        string toCompressedForm();

	bool operator==(const Coordinate& eqTo) {
		return X == eqTo.X && Y == eqTo.Y;
	}

    bool isInfinite() {
    	return X == 0 && Y == 0;
    }

	mpz_class X, Y;
};

ostream& operator<<(ostream& out, Coordinate& rhs);

class Jacobian {
public:

	Jacobian() {}
	Jacobian(mpz_class _x, mpz_class _y, mpz_class _z):
		X(_x), Y(_y), Z(_z) {}
	Jacobian(const Coordinate& rhs): //TODO: is this really the correct conversion?
                X(rhs.X), Y(rhs.Y), Z(1) {}

	static Jacobian infinity(){return Jacobian(1, 1, 0);}
    bool isInfinite() {
    	return X == 1 && Y == 1 && Z == 0;
    }
	mpz_class X, Y, Z;

};

/**
 * Lopez-Dahab projective coordinate
 * (X;Y;Z) corresponds to (X/Z;X/Z^2) in standard
 * coordinates
 *
 * See p.93, Menezes et. al.
 */
class LD {
public:

	LD() {}
	LD(mpz_class _x, mpz_class _y, mpz_class _z):
		X(_x), Y(_y), Z(_z) {}
	LD(const Coordinate& rhs): //TODO: is this really the correct conversion?
                X(rhs.X), Y(rhs.Y), Z(1) {}

	static LD infinity(){return LD(1, 0, 0);}
    bool isInfinite() {
    	return X == 1 && Y == 0 && Z == 0;
    }

	mpz_class X, Y, Z;

};

#endif /* COORDINATES_H_ */
