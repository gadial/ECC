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
	Ellipticcurve(mpz_class _mod, mpz_class _ECC_a, mpz_class _ECC_b):
		mod(_mod), ECC_a(_ECC_a), ECC_b(_ECC_b) {};
	Ellipticcurve(const char* _mod, int _mod_base,
			const char* _order, int _order_base,
			const char* _ecc_a, int _ecc_a_base,
			const char* _ecc_b, int _ecc_b_base,
			const char* _px, int _px_base,
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

	virtual mpz_class getOrder() {}

        //returns the point at infinity, as is represented by this class in the context of elliptic curves
        static Coordinate infinity(){return Coordinate(0,0);}

        static Ellipticcurve randomCurve(int number_of_bits, RandomNumberGenerator gen);
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
	 * Subtraction P-Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	Jacobian subtraction(Jacobian P, Coordinate Q);

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

private:

	/**
	 * Returns the non-adjacent form (NAF)
	 * of a positive integer k
	 */
	std::vector<int> getNAF(mpz_class k);

};

#endif /* ELLIPTICCURVE_H_ */
