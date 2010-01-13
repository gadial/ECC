/*
 * ECPrime.h
 *
 *  Created on: Nov 22, 2009
 *      Author: bhess
 */

#ifndef ECPRIME_H_
#define ECPRIME_H_

#include "ellipticcurve.h"
#include "zp_int.h"
class ECPrime: public Ellipticcurve {
public:
	ECPrime();
	ECPrime(mpz_class _mod, mpz_class _ECC_a, mpz_class _ECC_b) :
		Ellipticcurve(_mod, _ECC_a, _ECC_b) {}
	ECPrime(const char* _mod, int _mod_base, const char* _order,
			int _order_base, const char* _ecc_a, int _ecc_a_base,
			const char* _ecc_b, int _ecc_b_base, const char* _px, int _px_base,
			const char* _py, int _py_base) :
		Ellipticcurve(_mod, _mod_base, _order, _order_base, _ecc_a, _ecc_a_base,
			_ecc_b, _ecc_b_base, _px, _px_base, _py,  _py_base) {}


	static ECPrime randomCurve(int number_of_bits,
			RandomNumberGenerator gen);

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

	virtual ~ECPrime();

private:
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

};

#endif /* ECPRIME_H_ */
