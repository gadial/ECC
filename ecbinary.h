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
#include <iostream>

class ECBinary: public Ellipticcurve {
public:
	ECBinary();
	ECBinary(mpz_class _mod_poly, mpz_class _ECC_a, mpz_class _ECC_b) :
		Ellipticcurve(0, _ECC_a, _ECC_b) {

		// mod=2^_mod_2b
		//mod = 1;
		//mod <<= _mod_2b.get_ui();
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

private:
	LD doubling(LD P);
	LD addition(LD P, Coordinate Q);
	LD subtraction(LD P, Coordinate Q);

	mpz_class binMult(mpz_class a, mpz_class b, mpz_class f) {
		mpz_class res = 0;
		mpz_class const_1 = 2;
		mpz_class hiBitSet;
		for (int i = 0; i < mod; ++i) {
			if (b % const_1 == 1) {
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
};

#define WORD 32

/*
 * Implementation of the galois field
 * operations
 */
class GFE {
public:
	GFE(mpz_class _element, const mpz_class& _mod): mod(_mod), element(_element) {}

	/*
	 * Addition is just XOR
	 */
	GFE operator+(const GFE& other) {
		mpz_class res;
		mpz_xor(res.get_mpz_t(), element.get_mpz_t(), other.element.get_mpz_t());
		return GFE(res, mod);
	}

	/*
	 * Subtraction is the same as Addition
	 */
	GFE operator+(const GFE& other) {
		mpz_class res;
		mpz_xor(res.get_mpz_t(), element.get_mpz_t(), other.element.get_mpz_t());
		return GFE(res, mod);
	}

	/**
	 * Implements shift-and-add, complexity O(n) with degree n polynomial
	 * According to Lopez-Deneb (High-speed software mult. in F2m)
	 */
	GFE operator*(const GFE& other) {
		// TODO: more efficient impl. or use gf2x http://wwwmaths.anu.edu.au/~brent/gf2x.html
		mpz_class c = 0;
		int m = mpz_sizeinbase(mod.get_mpz_t(), 2);
		mpz_class mask_1 = 1 << (m - 1);
		mpz_class mask_2 = mask_1;
		mpz_class tmp;
		//int s = (m % WORD == 0 ? m / WORD : m / WORD + 1);
		//int k = m - 1 - WORD*(s - 1);
		for (int i = m - 1; i >= 0; --i, mask_1 >>= 1) {
			c <<= 1;
			mpz_and(tmp.get_mpz_t(), element.get_mpz_t(), mask_1.get_mpz_t());
			if (tmp != 0) {
				c ^= other.element;
			}
			mpz_and(tmp.get_mpz_t(), c.get_mpz_t(), mask_2.get_mpz_t());
			if (tmp != 0) {
				c ^= mod;
			}
		}
		return GFE(c, mod);

	}

	/*
	 * Returns the inverse Inverse
	 * Uses the Extended Euxlidean Algorithm over finite fields
	 */
	GFE operator!() {
		mpz_class u = element, v = mod;
		mpz_class g1 = 1, g2 = 0;
		while (u != 1) {
			int j = mpz_sizeinbase(u.get_mpz_t(), 2) - mpz_sizeinbase(v.get_mpz_t(), 2);
			if (j < 0) {
				mpz_class swap = u;
				u = v;
				v = swap;
				swap = g1;
				g1 = g2;
				g2 = swap;
				j = -j;
			}
			mpz_class sh = v << j;
			u ^= (v << j);
			g1 ^= (g2 << j);
		}
		return GFE(g1, mod);
	}

	void print() {
		std::cout << element.get_str(2) << " mod " << mod.get_str(2) << std::endl;
	}

	mpz_class mod;
	mpz_class element;
};

#endif /* ECBINARY_H_ */
