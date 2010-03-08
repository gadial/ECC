/*
 * ModPoly.cpp
 *
 *  Created on: Jan 27, 2010
 *      Author: bhess
 */

#include "ModPoly.h"

ModPoly::ModPoly(mpz_class bin) {
	int size = mpz_sizeinbase(bin.get_mpz_t(), 2);
	poly.rep.SetLength(size);
	mpz_class c1 = 1;
	for (int i = 0; i < size; ++i) {
		mpz_class sh = (bin >> i) & c1;
		poly.rep[i] = (sh == 1 ? 1 : 0);
	}
}

void ModPoly::set_precision(int N) {
	ZZ_p::init(to_ZZ(1) << N);
}

void ModPoly::set_mod(int deg, vector<int> pos, vector<int> val) {
	set_precision(deg);
	ZZ_pX P;
	P.rep.SetLength(deg + 1);
	for (unsigned int i = 0; i < pos.size(); ++i) {
		//SetCoeff(P, pos[i], val[i]);
		P.rep[pos[i]] = val[i];
	}
	ZZ_pE::init(P);
	cout << "modulus: " << ZZ_pE::modulus() << endl;
}

void ModPoly::set_coeff(int i, ZZ_p val) {
	poly.rep[i] = val;
}

ModPoly::ModPoly(int deg) {
	poly.rep.SetLength(deg + 1);
}

ModPoly ModPoly::operator+(const ModPoly& other) {
	return ModPoly(poly + other.poly);
}

ModPoly ModPoly::operator-(const ModPoly& other) {
	return ModPoly(poly - other.poly);
}

ModPoly ModPoly::operator*(const ModPoly& other) {
	return ModPoly(poly * other.poly);
}

ModPoly ModPoly::operator/(const ModPoly& other) {
	return ModPoly(poly / other.poly);
}

void ModPoly::operator *=(int s) {
	poly *= s;
}

void ModPoly::operator /=(int s) {

	for (int i = 0; i < poly.rep.length(); ++i) {
		poly.rep[i].LoopHole().rep[0].LoopHole() /= s;
		//cout << "representation: " << poly.rep[i].LoopHole().rep[0].LoopHole() << endl;
	}
	//poly /= s;
}

ModPoly ModPoly::operator <<(int m) {
	return ModPoly(poly << m);
}

ModPoly ModPoly::operator >>(int m) {
	return ModPoly(poly >> m);
}

void ModPoly::operator <<=(int m) {
	poly <<= m;
}

void ModPoly::operator >>=(int m) {
	poly >>= m;
}

bool ModPoly::operator==(const ModPoly& other) {
	return poly == other.poly;
}

ModPoly ModPoly::reverse() {
	return ModPoly(NTL::reverse(poly));
}

mpz_class ModPoly::to_gfe_el() {
	mpz_class res = 0;
	for (int i = 0; i < poly.rep.length(); ++i) {
		if (poly.rep[i].degree() % 2 == 1) {
			res |= (1 << i);
		}
	}
	return res;
}

ModPoly::~ModPoly() {
	// TODO Auto-generated destructor stub
}
