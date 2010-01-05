/*
 * Poly.cpp
 *
 *  Created on: Dec 19, 2009
 *      Author: bhess
 */

#include "Poly.h"

Poly::Poly(int _d) {
	coeffs = vector<mpz_class>(_d + 1, 0);
	degree = _d;
}

Poly::Poly(mpz_class bin) {
	int size = mpz_sizeinbase(bin.get_mpz_t(), 2);
	coeffs = vector<mpz_class>(size);
	degree = size - 1;
	for (int i = 0; i < size; ++i) {
		coeffs[i] = (bin >> i) & 1;
	}
}

void Poly::set_coeff(int d, mpz_class val) {
	coeffs[d] = val;
}

Poly Poly::operator+(const Poly& other) {
	Poly res = Poly(max(degree, other.degree));
	for (int i = 0; i < degree + 1 || i < other.degree + 1; ++i) {
		res.coeffs[i] = (i <= degree ? coeffs[i] : 0) + (i <= other.degree ? other.coeffs[i] : 0);
	}
	return res;
}

Poly Poly::operator-(const Poly& other) {
	Poly res = Poly(max(degree, other.degree));
	for (int i = 0; i < degree + 1 || i < other.degree + 1; ++i) {
		res.coeffs[i] = (i <= degree ? coeffs[i] : 0) - (i <= other.degree ? other.coeffs[i] : 0);
	}
	return res;
}

void Poly::operator*=(int s) {
	for (int i = 0; i < degree + 1; ++i) {
		coeffs[i] = coeffs[i] * s;
	}
}

void Poly::operator/=(int s) {
	for (int i = 0; i < degree + 1; ++i) {
		coeffs[i] = coeffs[i] / s;
	}
}

Poly Poly::operator*(const Poly& other) {
	Poly res = Poly(degree + other.degree);
	for (int i = 0; i < degree + 1; ++i) {
		for (int j = 0; j < other.degree + 1; ++j) {
			res.coeffs[i+j] += (coeffs[i] * other.coeffs[j]);
		}
	}
	return res;
}

void Poly::operator%=(mpz_class m) {
	for (int i = 0; i < degree + 1; ++i) {
		mpz_mod(coeffs[i].get_mpz_t(), coeffs[i].get_mpz_t(), m.get_mpz_t());
		//int aa = coeffs[i] % m;
		//coeffs[i] = (aa < 0 ? m + aa : aa);
	}
}

Poly Poly::operator>>(int m) {
	Poly res = Poly(degree - m);
	for (int i = m; i < degree + 1; ++i) {
		res.set_coeff(i - m, coeffs[i]);
	}
	return res;
}

Poly Poly::operator<<(int m) {
	Poly res = Poly(degree + m);
	for (int i = 0; i < degree + 1; ++i) {
		res.set_coeff(i + m, coeffs[i]);
	}
	return res;
}

Poly Poly::operator-() {
	Poly res(0);
	return res - (*this);
}

Poly Poly::PXpPmX() {
	Poly res = Poly(degree);
	for (int i = 0; i < degree + 1; i += 2) {
		res.coeffs[i] = coeffs[i] * 2;
	}
	for (int i = 1; i < degree + 1; i += 2) {
		res.coeffs[i] = 0;
	}
	return res;
}

bool Poly::operator==(const Poly& other) {
	for (int i = 0; i < (int)coeffs.size() && i < (int)other.coeffs.size(); ++i) {
		if (coeffs[i] != other.coeffs[i]) {
			return false;
		}
	}
	if (coeffs.size() > other.coeffs.size()) {
		for (int i = (int)other.coeffs.size(); i < (int)coeffs.size(); ++i) {
			if (coeffs[i] != 0) {
				return false;
			}
		}
	} else if (other.coeffs.size() > coeffs.size()) {
		for (int i = (int)coeffs.size(); i < (int)other.coeffs.size(); ++i) {
			if (other.coeffs[i] != 0) {
				return false;
			}
		}
	}
	return true;
}

Poly Poly::PXmPmX() {
	Poly res = Poly(degree);
	res.coeffs[0] = coeffs[0] * 2;
	for (int i = 1; i < degree + 1; i += 2) {
		res.coeffs[i] = coeffs[i] * 2;
	}
	for (int i = 2; i < degree + 1; i += 2) {
		res.coeffs[i] = 0;
	}
	return res;
}

Poly Poly::sqSubst() {
	Poly res = Poly(degree / 2);
	res.set_coeff(0, coeffs[0]);
	for (int i = 1; 2*i <= degree; ++i) {
		res.set_coeff(i, coeffs[2 * i]);
	}
	return res;
}

Poly Poly::reverse() {
	Poly res = Poly(degree);
	for (int i = 0; i < degree + 1; ++i) {
		res.set_coeff(i, coeffs[degree - i]);
	}
	return res;
}

Poly Poly::modXpowm(int m) {
	Poly res(m - 1);
	for (int i = 0; i < res.degree + 1; ++i) {
		res.set_coeff(i, coeffs[i]);
	}
	return res;
}

mpz_class Poly::to_gfe_el() {
	mpz_class res = 0;
	for (int i = 0; i < degree + 1; ++i) {
		if (coeffs[i] % 2 == 1) {
			res |= (1 << i);
		}
	}
	return res;
}

void Poly::print() {
	//cout << degree << endl;
	for (int i = degree; i >= 0; --i) {
		cout << coeffs[i] << ",";
	}
	cout << endl;
}
