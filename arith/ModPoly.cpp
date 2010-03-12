/*
 * ModPoly.cpp
 *
 *  Created on: Jan 27, 2010
 *      Author: bhess
 */

#include "ModPoly.h"

ModPoly::ModPoly(mpz_class bin, ZZ_pX _mod) {
	mod = _mod;
	int size = mpz_sizeinbase(bin.get_mpz_t(), 2);
	poly.rep.SetLength(size);
	mpz_class c1 = 1;
	for (int i = 0; i < size; ++i) {
		mpz_class sh = (bin >> i) & c1;
		//poly.rep[i] = (sh == 1 ? 1 : 0);
		SetCoeff(poly, i, to_ZZ_p((sh == 1 ? 1 : 0)));
	}
}

void ModPoly::set_precision(int N) {
	ZZ_p::init(to_ZZ(1) << N);
	ZZ_pX newPoly;
	ZZ_pX newMod;
	for (int i = 0; i < poly.rep.length(); ++i) {
		SetCoeff(newPoly, i, to_ZZ_p(poly.rep[i].LoopHole()));
	}
	for (int i = 0; i < mod.rep.length(); ++i) {
		SetCoeff(newMod, i, to_ZZ_p(mod.rep[i].LoopHole()));
	}
	// freeing current memory, otherwise NTL makes problems...
	kill();
	mod = newMod;
	poly = newPoly;
}

void ModPoly::set_mod(int deg, vector<int> pos, vector<int> val) {
	//set_precision(deg);
	mod.rep.SetLength(deg + 1);
	for (unsigned int i = 0; i < pos.size(); ++i) {
		//SetCoeff(P, pos[i], val[i]);
		SetCoeff(mod, pos[i], to_ZZ_p(val[i]));
		//P.rep[pos[i]] = val[i];
	}
	//ZZ_pE::init(P);
	//cout << "modulus: " << ZZ_pE::modulus() << endl;
}

void ModPoly::set_coeff(int i, ZZ_p val) {
	SetCoeff(poly, i, val);
}

ModPoly::ModPoly(int deg, ZZ_pX _mod) {
	mod = _mod;
	poly.rep.SetLength(deg + 1);
}

ModPoly ModPoly::operator+(const ModPoly& other) {
	//cout << "mod:" << mod << endl;
	//cout << poly << " + " << other.poly << " = ";
	ModPoly res = ModPoly((poly + other.poly), mod);
	if (res.poly != ZZ_pX::zero()) {
		res.poly = res.poly % mod;
	}
	//cout << res.poly << endl;
	return res;
}

ModPoly ModPoly::operator-(const ModPoly& other) {
	//cout << "mod:" << mod << endl;

	//cout << poly << " - " << other.poly << " = ";
	ModPoly res = ModPoly((poly - other.poly), mod);
	if (res.poly != ZZ_pX::zero()) {
		res.poly = res.poly % mod;
	}
	//cout << res.poly << endl;
	return res;
}

ModPoly ModPoly::operator*(const ModPoly& other) {
	//cout << poly << " * " << other.poly << " = ";
	ModPoly res = ModPoly((poly * other.poly), mod);
	if (res.poly != ZZ_pX::zero()) {
		res.poly = res.poly % mod;
	}
	//cout << res.poly << endl;
	return res;
}

ModPoly ModPoly::operator/(const ModPoly& other) {
	//cout << "mod:" << mod << endl;

	//cout << poly << " / " << other.poly << " = ";
	ModPoly res = ModPoly((poly / other.poly), mod);
	if (res.poly != ZZ_pX::zero()) {
		res.poly = res.poly % mod;
	}
	//cout << res.poly << endl;
	return res;
}

void ModPoly::operator *=(int s) {
	//cout << poly << " * " << s << " = ";
	poly = (poly * s);
	if (poly != ZZ_pX::zero())
		poly = poly % mod;
	//cout << poly << endl;
}

void ModPoly::operator /=(int s) {
	//cout << poly << " / " << s << " = ";
	for (int i = 0; i < poly.rep.length(); ++i) {
		SetCoeff(poly, i, to_ZZ_p(poly.rep[i].LoopHole() / s));
		//poly.rep[i].LoopHole().rep[0].LoopHole() /= s;
		//cout << "representation: " << poly.rep[i].LoopHole().rep[0].LoopHole() << endl;
	}
	//cout << poly << endl;
	//poly /= s;
}

ModPoly ModPoly::operator <<(int m) {
	ModPoly res = ModPoly((poly << m) % mod, mod);
	return res;
}

ModPoly ModPoly::operator >>(int m) {
	ModPoly res = ModPoly((poly >> m) % mod, mod);
	return res;
}

void ModPoly::operator <<=(int m) {
	poly = (poly << m) % mod;
}

void ModPoly::operator >>=(int m) {
	poly = (poly >> m) % mod;
}

bool ModPoly::operator==(const ModPoly& other) {
	return poly == other.poly;
}

ModPoly ModPoly::reverse() {
	ModPoly res = ModPoly(NTL::reverse(poly), mod);
	return res;
}

mpz_class ModPoly::to_gfe_el() {
	mpz_class res = 0;
	for (int i = 0; i < poly.rep.length(); ++i) {
		if (poly.rep[i].LoopHole() % to_ZZ(2) == 1) {
			res |= (1 << i);
		}
	}
	return res;
}

ModPoly::~ModPoly() {
	// TODO Auto-generated destructor stub
}
