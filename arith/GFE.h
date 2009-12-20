/*
 * GFE.h
 *
 *  Created on: Dec 19, 2009
 *      Author: bhess
 */

#ifndef GFE_H_
#define GFE_H_

#include <iostream>
#include <gmpxx.h>

/*
 * Implementation of the galois field
 * operations
 */
class GFE {
public:
	// Dummy constructor.
	GFE() {}

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
	GFE operator-(const GFE& other) {
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

	/**
	 * degree of the modulus
	 */
	int mod_deg() {
		return mpz_sizeinbase(mod.get_mpz_t(), 2) - 1;
	}

	/*
	 * degree of the element
	 */
	int el_deg() {
		return mpz_sizeinbase(element.get_mpz_t(), 2) - 1;
	}

	/**
	 * is the element monic
	 * i.e. element has highes possible degree
	 */
	bool is_monic() {
		return mod_deg() == el_deg();
	}

	mpz_class mod;
	mpz_class element;
};

#endif /* GFE_H_ */
