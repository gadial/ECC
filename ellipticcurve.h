/*
 * ellipticcurve.h
 *
 *  Created on: Nov 5, 2009
 *      Author: bhess
 */

#ifndef ELLIPTICCURVE_H_
#define ELLIPTICCURVE_H_

#include <vector>
#include <gmpxx.h>
#include "coordinates.h"

class Ellipticcurve {
public:
	Ellipticcurve();
	virtual ~Ellipticcurve();

	/*
	 * --------------------------------
	 * Definition of the Elliptic Curve
	 * --------------------------------
	 *
	 * The prime modulus for the underlying finite field
	 */
	mpz_class mod;

	/**
	 * Equation that defines the EC:
	 * y^2 = x^3 + ax + b
	 */
	int ECC_a, ECC_b;

	/**
	 * A point in E(F_p)
	 */
	Coordinate point;

	/*
	 * -------------
	 * EC operations
	 * -------------
	 */

	/**
	 * Addition P+Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	Jacobian addition(Jacobian P, Coordinate Q);

	/**
	 * Doubling of a point P -> 2P
	 * in jacobian coordinates
	 */
	Jacobian doubling(Jacobian P);

	/**
	 * Repeated doubling of a point P in jacobian coordinates
	 * (m times) -> 2^m P
	 */
	Jacobian repeatedDoubling(Jacobian P, int m);

	/**
	 * Point multiplication (k times)
	 * -> kP
	 */
	Jacobian pointMultiplication(Coordinate P, mpz_class k);

private:

	/**
	 * Returns the non-adjacent form (NAF)
	 * of a positive integer k
	 */
	std::vector<int> getNAF(mpz_class k);

	/**
	 * Returns the point -P
	 * According to p.80 (char /= 2,3)
	 */
	Coordinate getNegative(const Coordinate& P);
};

#endif /* ELLIPTICCURVE_H_ */
