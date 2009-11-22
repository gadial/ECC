/*
 * ECBinary.h
 *
 *  Created on: Nov 22, 2009
 *      Author: bhess
 *
 * Representation of an elliptic curve over a finite
 * field with characteristic 2: F_2^m
 *
 * Equation: (non-supersingular)
 * y^2+xy=x^3+Ax^2+b
 *
 */

#ifndef ECBINARY_H_
#define ECBINARY_H_

#include "ellipticcurve.h"

class ECBinary: public Ellipticcurve {
public:
	ECBinary();
	ECBinary(mpz_class _mod_2b, mpz_class _ECC_a, mpz_class _ECC_b) :
		Ellipticcurve(0, _ECC_a, _ECC_b) {

		// mod=2^_mod_2b
		mod = 1;
		mod <<= _mod_2b.get_ui();
	}
	ECBinary(const char* _mod_2b, int _mod_base, const char* _order,
			int _order_base, const char* _ecc_a, int _ecc_a_base,
			const char* _ecc_b, int _ecc_b_base, const char* _px, int _px_base,
			const char* _py, int _py_base) :
		Ellipticcurve(_mod_2b, _mod_base, _order, _order_base, _ecc_a,
				_ecc_a_base, _ecc_b, _ecc_b_base, _px, _px_base, _py, _py_base) {

		// mod = 2^_mod_2b
		unsigned long int tmp = mod.get_ui();
		mod = 1;
		mod <<= tmp;
	}

	/**
	 * Addition P+Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	Coordinate addition(Coordinate P, Coordinate Q);

	/**
	 * Subtraction P-Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	Coordinate subtraction(Coordinate P, Coordinate Q);

	/**
	 * Doubling of a point P -> 2P
	 * in jacobian coordinates
	 */
	Coordinate doubling(Coordinate P);

	/**
	 * Repeated doubling of a point P in jacobian coordinates
	 * (m times) -> 2^m P
	 */
	Coordinate repeatedDoubling(Coordinate P, int m);

	/**
	 * Point multiplication (k times)
	 * -> kP
	 */
	Coordinate pointMultiplication(Coordinate P, mpz_class k);

private:
	LD doubling(LD P);
	LD addition(LD P, Coordinate Q);
	LD subtraction(LD P, Coordinate Q);
};

#endif /* ECBINARY_H_ */
