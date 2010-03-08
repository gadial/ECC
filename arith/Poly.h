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

	static Poly zero() {
		Poly a(0);
		a.set_coeff(0, 0);
		return a;
	}

	static Poly one() {
		Poly a(0);
		a.set_coeff(0, 1);
		return a;
	}

	static Poly two() {
		Poly a(0);
		a.set_coeff(0, 2);
		return a;
	}

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
	//void operator%=(mpz_class m);
	void operator%=(int logm);
	Poly operator>>(int m);
	Poly operator<<(int m);
	Poly operator-();
	bool operator==(const Poly& other);

	Poly PXpPmX();
	Poly PXmPmX();
	Poly sqSubst();
	Poly modXpowm(int m);

	Poly reverse();

	mpz_class to_gfe_el();

	void print();
};
#endif /* POLY_H_ */
