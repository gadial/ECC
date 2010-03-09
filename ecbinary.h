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
 * y^2+xy=x^3+Ax^2+B
 *
 */

#ifndef ECBINARY_H_
#define ECBINARY_H_

#include "ellipticcurve.h"
#include "arith/GFE.h"
//#include <iostream>


//class GFE;

class ECBinary: public Ellipticcurve {
public:
	ECBinary();
	ECBinary(mpz_class _mod_poly, mpz_class _ECC_a, mpz_class _ECC_b) :
		Ellipticcurve(_mod_poly, _ECC_a, _ECC_b) {

		// mod=2^_mod_2b
		//mod = 1;
		//mod <<= _mod_2b.get_ui();
	}
	ECBinary(int _tri_d, int _tri_k, mpz_class _ECC_a, mpz_class _ECC_b) :
			Ellipticcurve(0, _ECC_a, _ECC_b) {

		tri_d = _tri_d;
		tri_k = _tri_k;
		mod |= (1 << tri_d);
		mod |= (1 << tri_k);
		mod |= 1;

		sqrtx = GFE::get_sqrtx(tri_d, tri_k, mod);


		sqrtx.print();
		GFE xx = sqrtx * sqrtx;
		xx.print();

		GFE el(0b101011, mod);
		el.print();
		GFE elsq = el * el;
		elsq.print();
		GFE back = elsq.get_sqrt(sqrtx);
		back.print();

	}
	ECBinary(const char* _mod_poly, int _mod_base, const char* _order,
			int _order_base, const char* _ecc_a, int _ecc_a_base,
			const char* _ecc_b, int _ecc_b_base, const char* _px, int _px_base,
			const char* _py, int _py_base) :
		Ellipticcurve(_mod_poly, _mod_base, _order, _order_base, _ecc_a,
				_ecc_a_base, _ecc_b, _ecc_b_base, _px, _px_base, _py, _py_base) {

		// mod = 2^_mod_2b
		//unsigned long int tmp = mod.get_ui();
		//mod = 1;
		//mod <<= tmp;
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

	std::vector<GFE> get_jinvariants();

	std::vector<GFE> update_js(int n, std::vector<GFE> Jinvs);

	std::vector<GFE> lift_jinvariants(int n, std::vector<GFE> jinvs);

	mpz_class satohfgh_point_counting();

	Coordinate getPointCompressedForm(string from) {
		// TODO: implement...
	}

	/*
	 * Gets a point, given the x coordinate for the following equation:
	 *  y^2 + xy = x^3 + ax^2 + b
	 * Solves the quadratic equation f(y) = y^2 + xy - (x^3 + ax^2 + b) = 0
	 */
	Coordinate getPoint_interface(mpz_class x, bool negative_value = false);

private:

	/*
	 * modular equation phi_2
	 */
	inline GFE phi_2(GFE x, GFE y);

	/*
	 * derivation of phi_2 with respect to x
	 */
	inline GFE phi_2_x(GFE x, GFE y);

	inline mpz_class pi_j(mpz_class J) {
		return 0;
	}

	LD doubling(LD P);
	LD addition(LD P, Coordinate Q);
	LD subtraction(LD P, Coordinate Q);

	mpz_class binMult(mpz_class a, mpz_class b, mpz_class f) {
		mpz_class res = 0;
		mpz_class const_2 = 2;
		mpz_class hiBitSet;
		for (int i = 0; i < mod; ++i) {
			if (b % const_2 == 1) {
				res ^= a;
			}
			mpz_and(hiBitSet.get_mpz_t(), a.get_mpz_t(), mod.get_mpz_t());
			a <<= 1;
			if (hiBitSet != 0) {
				a ^= f;
			}
			b >>= 1;
		}
		return res;
	}

	Coordinate toCoordinate(const LD& ld);

	int tri_d, tri_k;
	GFE sqrtx;
};

#define WORD 32



#endif /* ECBINARY_H_ */
