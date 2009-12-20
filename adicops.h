/*
 * adicops.h
 *
 *  Created on: Dec 15, 2009
 *      Author: bhess
 */

#ifndef ADICOPS_H_
#define ADICOPS_H_

#include <vector>
#include <gmpxx.h>
#include "arith/Poly.h"
using namespace std;


class adicops {
public:
	adicops();
	void do_s();

	/**
	 * Computes the Teichmüller modulus
	 * with precision N
	 */
	Poly get_teichmuller_modulus(Poly in, int N);
	Poly teichmuller_modulus_increment(const Poly& M0, const Poly& M1,
			const Poly& V, int N);
	Poly poly_remainder(Poly a, Poly b, int prec);
	Poly poly_invert(Poly f, int p, int prec);
	Poly get_inverse(Poly a, Poly mod, int prec);
	Poly get_invsqrt(Poly a, Poly approx, Poly mod, int prec);
	Poly get_sqrt(Poly a, Poly approx, Poly mod, int prec);


	mpz_class get_points(mpz_class _c, mpz_class _mod, int d);

	virtual ~adicops();
};



#endif /* ADICOPS_H_ */
