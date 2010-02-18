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
#include "primes.h"

class Ellipticcurve {
public:
	Ellipticcurve();
	Ellipticcurve(mpz_class _mod, mpz_class _ECC_a, mpz_class _ECC_b) :
		mod(_mod), ECC_a(_ECC_a), ECC_b(_ECC_b) {
	}
	;
	Ellipticcurve(const char* _mod, int _mod_base, const char* _order,
			int _order_base, const char* _ecc_a, int _ecc_a_base,
			const char* _ecc_b, int _ecc_b_base, const char* _px, int _px_base,
			const char* _py, int _py_base);

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
	mpz_class ECC_a, ECC_b;

	/**
	 * A point in E(F_p)
	 */
	Coordinate point;

	int get_bits();

	void set_point_compressed(string cpoint) {
		point = getPointCompressedForm(cpoint);
	}

	mpz_class getOrder() {
		return order;
	}

    void setOrder(mpz_class _order) {
            order = _order;
	}

	//static Ellipticcurve randomCurve(int number_of_bits,
	//		RandomNumberGenerator gen);

	/*
	 * -------------
	 * EC operations
	 * -------------
	 */

	/**
	 * Addition P+Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	virtual Coordinate addition(Coordinate P, Coordinate Q) = 0;

	/**
	 * Subtraction P-Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	virtual Coordinate subtraction(Coordinate P, Coordinate Q) = 0;

	/**
	 * Doubling of a point P -> 2P
	 * in jacobian coordinates
	 */
	virtual Coordinate doubling(Coordinate P) = 0;

	/**
	 * Repeated doubling of a point P in jacobian coordinates
	 * (m times) -> 2^m P
	 */
	virtual Coordinate repeatedDoubling(Coordinate P, int m) = 0;

	/**
	 * Point multiplication (k times)
	 * -> kP
	 */
	virtual Coordinate pointMultiplication(Coordinate P, mpz_class k) = 0;

	/*
	 * -----------
	 * Accesseors
	 * ----------
	 *
	 * /

	 /**
	 *
	 * Given x, returns the point corresponding to (x,y) (or (x,-y) if asked)
	 * If there is no corresponding point, returns the infinity
	 */
	Coordinate getPoint(mpz_class x, bool negative_value = false);

	/*
	 * Gets the point from compressed format
	 * see IEEE P1363 / D8 E.2.3.1
	 */
	virtual Coordinate getPointCompressedForm(string from) = 0;
protected:

	/*
	 * Order of the EC
	 */
	mpz_class order;

	/**
	 * Returns the point -P
	 * According to p.80 (char /= 2,3)
	 */
	Coordinate getNegative(const Coordinate& P);

	/**
	 * Returns the non-adjacent form (NAF)
	 * of a positive integer k
	 */
	std::vector<int> getNAF(mpz_class k);

};

#endif /* ELLIPTICCURVE_H_ */
