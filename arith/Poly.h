/*
 * Poly.h
 *
 *  Created on: Dec 19, 2009
 *      Author: bhess
 */

#ifndef POLY_H_
#define POLY_H_

#include <vector>
#include <gmpxx.h>
#include "GFE.h"

using namespace std;

class Poly {
public:
	Poly() {};
	Poly(mpz_class bin);
	Poly(int d);

	void set_coeff(int d, mpz_class val);

	vector<mpz_class> coeffs;
	int degree;

	Poly operator+(const Poly& other);
	Poly operator-(const Poly& other);
	Poly operator*(const Poly& other);
	void operator*=(int s);
	void operator/=(int s);
	void operator%=(mpz_class m);
	Poly operator>>(int m);
	Poly operator<<(int m);
	Poly operator-();

	Poly PXpPmX();
	Poly PXmPmX();
	Poly sqSubst();
	Poly modXpowm(int m);

	Poly reverse();

	mpz_class to_gfe_el();

	void print();
};
#endif /* POLY_H_ */
