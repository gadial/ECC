/*
 * GFE.cpp
 *
 *  Created on: Jan 24, 2010
 *      Author: bhess
 */

#include "GFE.h"

#ifdef USE_NTL

GFE::GFE(mpz_class _element, const mpz_class& _mod) {
	GF2X gx;
	convert(gx, _mod);
	init(gx);
	convert(ntl_element, _element);
}

GFE GFE::operator+(const GFE& other) {
	//std::cout << ntl_element << "+" << other.ntl_element << "=";
	GFE ret = GFE(ntl_element + other.ntl_element);
	//std::cout << ret.ntl_element << std::endl;
	return ret;
}

GFE GFE::operator-(const GFE& other) {
	return GFE(ntl_element - other.ntl_element);
}

GFE GFE::operator*(const GFE& other) {
	//std::cout << ntl_element << "*" << other.ntl_element << "=";
	GFE ret = GFE(ntl_element * other.ntl_element);
	//std::cout << ret.ntl_element << std::endl;
	return ret;
}

GFE GFE::operator!() {
	//std::cout << "!" << ntl_element << "=";
	GFE ret = GFE(inv(ntl_element));
	//std::cout << ret.ntl_element << std::endl;
	return ret;
}

void GFE::print() {
	std::cout << ntl_element << ", mod " << GF2E::modulus() << std::endl;
	//std::cout << element.get_str(2) << " mod " << mod.get_str(2) << std::endl;
}

/**
 * degree of the modulus
 */
int GFE::mod_deg() {
	return GF2E::degree();
}

/*
 * degree of the element
 */
int GFE::el_deg() {
	return deg(ntl_element.LoopHole());
}

/**
 * is the element monic
 * i.e. element has highes possible degree
 */
bool GFE::is_monic() {
	return mod_deg() == el_deg();
}

GFE GFE::get_sqrt() {

	//char* s = mod.get_;
	//const char* s = (const char*)mod.get_str(10);
	//NTL::ZZ mo = NTL::to_ZZ(s);


}

GFE GFE::get_sqrt(GFE sqrtx) {
	//GFE x(0b10, mod);
	/*
	 GFE even = GFE(element & 0b1, mod);
	 GFE odd = GFE((element & 0b10) >> 1, mod);
	 */

	vec_GF2 even; even.SetLength(el_deg() + 1);
	vec_GF2 odd; odd.SetLength(el_deg() + 1);

	vec_GF2 ell = to_vec_GF2(ntl_element.LoopHole());
	//std::cout << "ntl: " << ell << ", deg: " << el_deg() << std::endl;

	for (int i = 0; i <= el_deg() / 2; ++i) {
		if (ell[2 * i] == 1) {
			even[i] = 1;// |= (1 << i);
		}
		if (ell.length() > (2 * i + 1)) {
			if (ell[(2 * i) + 1] == 1) {
				odd[i] = 1;// |= (1 << i);
			}
		}

		/*
		if ((element >> (2 * i)) % 2 == 1) {
			even[i] = 1;// |= (1 << i);
		}
		if ((element >> ((2 * i) + 1)) % 2 == 1) {
			odd[i] = 1;// |= (1 << i);
		}
		*/
	}


	return GFE(even) + (sqrtx * GFE(odd));
}

/*
 * Solves T^2+aT+b=0, to T
 */
GFE GFE::solve_quad_eq(GFE a, GFE b) {
	return solve_quad_eq(b * !(a * a));
}

/*
 * Solves T^2 + T = c, to T
 */
GFE GFE::solve_quad_eq(GFE c) {

	// TODO: only for odd degree polynomials yet.
	int d = c.el_deg();
	GFE csq = c * c;
	GFE res = csq;
	// sum c^2^(2i+1), i=0..(d-3)/2
	// c^2^1+c^2^3+c^2^5+...=c^2+c^8+c^32+c^128
	// c^j=c^(j-1)^4=((c^(j-1))^2)^2
	for (int i = 1; i <= (d - 3) / 2; ++i) {
		// squaring two times, then add...
		csq = csq * csq;
		csq = csq * csq;
		res = res + csq;
	}
	return res;
}

GFE GFE::get_sqrtx(int d, int k, mpz_class mod) {
	GFE x(0b10, mod);
	GFE sqx_k = x;
	if (k % 2 == 1) {
		int i = 1;
		for (; i < (k + 1) / 2; ++i) {
			sqx_k = sqx_k * x;
		}
		GFE sqx_d = sqx_k;
		for (; i < (d + 1) / 2; ++i) {
			sqx_d = sqx_d * x;
		}

		return (sqx_k + sqx_d);
	} else {
		int i = 1;
		for (; i < k / 2; ++i) {
			sqx_k = sqx_k * x;
		}
		GFE sqx_d = sqx_k;
		for (; i < (d - 1) / 2; ++i) {
			sqx_d = sqx_d * x;
		}

		return ((!sqx_d) * (sqx_k + GFE(1, mod)));
	}
}

mpz_class GFE::get_element() {
	mpz_class res;
	convert(res, ntl_element);
	return res;
}

NTL::GF2E GFE::sqr(NTL::GF2E x, NTL::GF2X mod) {
	NTL::GF2E::init(mod);

}

#else

GFE::GFE(mpz_class _element, const mpz_class& _mod) {
	element = _element;
	mod = _mod;
}

GFE GFE::get_sqrtx(int d, int k, mpz_class mod) {
		GFE x(0b10, mod);
		GFE sqx_k = x;
		if (k % 2 == 1) {
			int i = 1;
			for (; i < (k + 1) / 2; ++i) {
				sqx_k = sqx_k * x;
			}
			GFE sqx_d = sqx_k;
			for (; i < (d + 1) / 2; ++i) {
				sqx_d = sqx_d * x;
			}

			return (sqx_k + sqx_d);
		} else {
			int i = 1;
			for (; i < k / 2; ++i) {
				sqx_k = sqx_k * x;
			}
			GFE sqx_d = sqx_k;
			for (; i < (d - 1) / 2; ++i) {
				sqx_d = sqx_d * x;
			}

			return ((!sqx_d) * (sqx_k + GFE(1, mod)));
		}
	}

	/*
	 * Addition is just XOR
	 */
	GFE GFE::operator+(const GFE& other) {
		mpz_class res;
		mpz_xor(res.get_mpz_t(), element.get_mpz_t(), other.element.get_mpz_t());
		return GFE(res, mod);
	}

	/*
	 * Subtraction is the same as Addition
	 */
	GFE GFE::operator-(const GFE& other) {
		mpz_class res;
		mpz_xor(res.get_mpz_t(), element.get_mpz_t(), other.element.get_mpz_t());
		return GFE(res, mod);
	}

	/**
	 * Implements shift-and-add, complexity O(n) with degree n polynomial
	 * According to Lopez-Deneb (High-speed software mult. in F2m)
	 */
	GFE GFE::operator*(const GFE& other) {
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
	GFE GFE::operator!() {
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

	GFE GFE::get_sqrt() {

		//char* s = mod.get_;
		//const char* s = (const char*)mod.get_str(10);
		//NTL::ZZ mo = NTL::to_ZZ(s);


	}

	GFE GFE::get_sqrt(GFE sqrtx) {
		//GFE x(0b10, mod);
		/*
		GFE even = GFE(element & 0b1, mod);
		GFE odd = GFE((element & 0b10) >> 1, mod);
		*/
		mpz_class even = 0;
		mpz_class odd = 0;

		for (int i = 0; i <= el_deg() / 2; ++i) {
			if ((element >> (2 * i)) % 2 == 1) {
				even |= (1 << i);
			}
			if ((element >> ((2 * i) + 1)) % 2 == 1) {
				odd |= (1 << i);
			}
		}

		return GFE(even, mod) + (sqrtx * GFE(odd, mod));
	}

	/*
	 * Solves T^2+aT+b=0, to T
	 */
	GFE GFE::solve_quad_eq(GFE a, GFE b) {
		return solve_quad_eq(b * !(a*a));
	}

	/*
	 * Solves T^2 + T = c, to T
	 */
	GFE GFE::solve_quad_eq(GFE c) {

		// TODO: only for odd degree polynomials yet.
		int d = c.el_deg();
		GFE csq = c * c;
		GFE res = csq;
		// sum c^2^(2i+1), i=0..(d-3)/2
		// c^2^1+c^2^3+c^2^5+...=c^2+c^8+c^32+c^128
		// c^j=c^(j-1)^4=((c^(j-1))^2)^2
		for (int i = 1; i <= (d - 3) / 2; ++i) {
			// squaring two times, then add...
			csq = csq * csq;
			csq = csq * csq;
			res = res + csq;
		}
		return res;
	}

	void GFE::print() {
		std::cout << element.get_str(2) << " mod " << mod.get_str(2) << std::endl;
	}

	/**
	 * degree of the modulus
	 */
	int GFE::mod_deg() {
		return mpz_sizeinbase(mod.get_mpz_t(), 2) - 1;
	}

	/*
	 * degree of the element
	 */
	int GFE::el_deg() {
		return mpz_sizeinbase(element.get_mpz_t(), 2) - 1;
	}

	/**
	 * is the element monic
	 * i.e. element has highes possible degree
	 */
	bool GFE::is_monic() {
		return mod_deg() == el_deg();
	}

	mpz_class GFE::get_element() {
		return element;
	}
#endif
