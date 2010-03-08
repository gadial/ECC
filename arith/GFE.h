/*
 * GFE.h
 *
 *  Created on: Dec 19, 2009
 *      Author: bhess
 */

#ifndef GFE_H_
#define GFE_H_

// Uncomment if you want to use NTL
//#define USE_NTL

#include <iostream>
#include <gmpxx.h>
#include <string>
#ifdef USE_NTL
#include <NTL/GF2E.h>
using namespace NTL;
#endif

/*
 * Implementation of the galois field
 * operations
 */
class GFE {
public:
	// Dummy constructor.
#ifdef USE_NTL
	static void convert(GF2E& to, const mpz_class& from) {
		GF2X x;
		convert(x, from);
		to = to_GF2E(x);
	}

	static void convert(GF2X& to, const mpz_class& from) {
		std::string sst = from.get_str(2);
		const char* cc = sst.c_str();
		vec_GF2 vgf2;
		vgf2.SetLength(sst.length());
		for (int i = 0; i < (int)sst.length(); ++i) {
			vgf2[i] = (int)cc[i];
		}
		to = to_GF2X(vgf2);
		//std::cout << vgf2 << std::endl;
	}

	static void convert(mpz_class& to, const GF2X& from) {
		vec_GF2 v = to_vec_GF2(from);
		mpz_class c1 = 1;
		to = 0;
		for (int i = 0; i < v.length(); ++i) {
			if (v[i] == 1) {
				to |= (c1 << i);
			}
		}
	}

	static void convert(mpz_class& to, GF2E& from) {
		convert(to, from.LoopHole());
	}

	static NTL::GF2E sqr(NTL::GF2E x, NTL::GF2X mod);

	/*
	 * Initialize with a modululus
	 */
	static void init(const GF2X& newMod) {
		GF2E::init(newMod);
	}
	static void init(const mpz_class& newMod) {
		std::cout << newMod << std::endl;
		GF2X g;
		convert(g, newMod);
		GF2E::init(g);
	}

	GFE(NTL::GF2E _element) { ntl_element = _element; }
	GFE(NTL::vec_GF2 _element) { ntl_element = to_GF2E(to_GF2X(_element)); }
	//GFE(NTL::ZZ _element) { conv(ntl_element, _element); }

	GF2E ntl_element;
#endif
	GFE() {}

	GFE(mpz_class _element, const mpz_class& _mod);

	/*
	 * solves x = a^2 to a
	 */
	static GFE get_sqrtx(int d, int k, mpz_class mod);

	/*
	 * Addition is just XOR
	 */
	GFE operator+(const GFE& other);

	/*
	 * Subtraction is the same as Addition
	 */
	GFE operator-(const GFE& other);

	/**
	 * Implements shift-and-add, complexity O(n) with degree n polynomial
	 * According to Lopez-Deneb (High-speed software mult. in F2m)
	 */
	GFE operator*(const GFE& other);

	bool operator==(const mpz_class& other) {
		return get_element() == other;
	}

	bool operator==(const GFE& other) {
		return get_element() == other.element;
	}

	/*
	 * Returns the inverse Inverse
	 * Uses the Extended Euxlidean Algorithm over finite fields
	 */
	GFE operator!();

	GFE get_sqrt(GFE sqrtx);

	GFE get_sqrt();

	/*
	 * Solves T^2+aT+b=0, to T
	 */
	static GFE solve_quad_eq(GFE a, GFE b);

	/*
	 * Solves T^2 + T = c, to T
	 */
	static GFE solve_quad_eq(GFE c);

	void print();

	/**
	 * degree of the modulus
	 */
	int mod_deg();

	/*
	 * degree of the element
	 */
	int el_deg();

	/**
	 * is the element monic
	 * i.e. element has highes possible degree
	 */
	bool is_monic();

	mpz_class mod;
	mpz_class get_element();
	mpz_class element;
};

#endif /* GFE_H_ */
