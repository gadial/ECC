/*
 * ModPoly.h
 *
 *  Created on: Jan 27, 2010
 *      Author: bhess
 */

#ifndef MODPOLY_H_
#define MODPOLY_H_

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>
#include <vector>
#include <gmpxx.h>

using namespace NTL;
using namespace std;

class ModPoly {
public:

	/*
	 * Set precision to 2^N
	 */
	static void set_precision(int N);

	/*
	 * set the modulus for the polynomial
	 * deg: degree
	 * pos: setting coefficients c_pos
	 * val: setting c_pos to val
	 */
	static void set_mod(int deg, vector<int> pos, vector<int> val);

	/*
	 * returns the polynomial 0
	 */
	static ModPoly zero() {
		return ModPoly(ZZ_pEX::zero());
	}

	/*
	 * returns the polynomial 1
	 */
	static ModPoly one() {
		return constP(1);
	}

	/*
	 * returns the constant polynomial i
	 */
	static ModPoly constP(int i) {
		ZZ_pEX p;
		p.rep.SetLength(1);
		p.rep[0] = i;
		return ModPoly(p);
	}

	ModPoly(int deg);
	ModPoly(mpz_class bin);
	ModPoly(ZZ_pEX _p): poly(_p) {}

	void set_coeff(int i, ZZ_p val);

	ModPoly operator+(const ModPoly& other);
	ModPoly operator-(const ModPoly& other);
	ModPoly operator*(const ModPoly& other);
	void operator*=(int s);
	void operator/=(int s);
	//void operator%=(mpz_class m);
	//void operator%=(int logm);
	ModPoly operator>>(int m);
	ModPoly operator<<(int m);
	void operator>>=(int m);
	void operator<<=(int m);
	ModPoly operator-();
	bool operator==(const ModPoly& other);
	ModPoly operator/(const ModPoly& other);

	ModPoly PXpPmX();
	ModPoly PXmPmX();
	ModPoly sqSubst();
	ModPoly modXpowm(int m);

	ModPoly reverse();

	mpz_class to_gfe_el();

	virtual ~ModPoly();

	ZZ_pEX poly;
};

#endif /* MODPOLY_H_ */
