/*
 * adicops.h
 *
 *  Created on: Dec 15, 2009
 *      Author: bhess
 */

#ifndef ADICOPS_H_
#define ADICOPS_H_

//#define VERBOSE

#include <vector>
#include <gmpxx.h>
#include "arith/Poly.h"
#include "arith/ModPoly.h"
using namespace std;

class Adicops {
public:
	Adicops(Poly mod);
	static void do_s();
	static void do_ntl();

	Poly get_mod();
	/**
	 * Computes the Teichmüller modulus
	 * with precision N
	 */
	void set_teichmuller_modulus(Poly in, int prec);
	Poly teichmuller_modulus_increment(const Poly& M0, const Poly& M1,
			const Poly& V, int N);
	Poly poly_remainder(Poly a, int prec);
	Poly poly_remainder(Poly a, Poly b, int prec);
	Poly poly_division_rem(Poly num, Poly denom, int prec);
	Poly poly_division(Poly num, Poly denom, int prec);
	Poly poly_invert(Poly f, int p, int prec);
	Poly get_inverse(Poly a, int prec);
	Poly get_invsqrt(Poly a, Poly approx, int prec);
	Poly get_sqrt(Poly a, int prec);

	ModPoly get_inverse(ModPoly a, int prec);
	ModPoly get_invsqrt(ModPoly a, ModPoly approx, int prec);
	ModPoly get_sqrt(ModPoly a, int prec);

	bool testsqrt(Poly sqrt, Poly in, Poly mod, int prec);


	mpz_class get_points_AGM_bivariate(mpz_class _c, int d);
	ZZ get_points_AGM_bivariate_v2(mpz_class _c, int d);
	mpz_class get_points_AGM_univariate(mpz_class _c, mpz_class _mod, int d);

private:
	Poly get_teichmuller_modulus(Poly in, int prec);
	Poly mod;
};



#endif /* ADICOPS_H_ */
