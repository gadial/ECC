/*
 * adicops.cpp
 *
 *  Created on: Dec 15, 2009
 *      Author: bhess
 */

#include "adicops.h"
#include <iostream>
#include <cassert>
#define NDEBUG

//#include <NTL/GF2X.h>

Adicops::Adicops(Poly _mod) {
	mod = _mod;
}

void Adicops::do_s() {

	int deg = 7;
	mpz_class c1 = 1;
	mpz_class orer = (c1 << deg);
	mpz_class tmp;
	tmp.set_str("3", 16);
	tmp |= orer;
	Adicops a(tmp);
	//GFE::init(tmp);
	vector<int> pos(3), val(3);
	pos[0] = 0; pos[1] = 6; pos[2] = 7;
	val[0] = 1; val[1] = 1; val[2] = 1;
	ModPoly::set_mod(deg, pos, val);
	cout << a.get_points_AGM_bivariate_v2(0b10001, deg) << " points!" << endl;
}

void Adicops::do_ntl() {

}

Poly Adicops::get_mod() {
	return mod;
}

void Adicops::set_teichmuller_modulus(Poly in, int prec) {
	mod = get_teichmuller_modulus(in, prec);
}

Poly Adicops::get_teichmuller_modulus(Poly in, int N) {
	Poly M;
	if (N == 1) {
		M = in;
		M %= 1;
#ifdef VERBOSE
		cout << "----------------" << endl << "N = " << N << endl;

		cout << "M: ";
		M.print();
#endif
	} else {
		int newN = (N % 2 == 0 ? N / 2 : N / 2 + 1);
		M = get_teichmuller_modulus(in, newN);
#ifdef VERBOSE
		cout << "----------------" << endl << "N = " << N << ", N' = " << newN << endl;
		cout << "M: ";

		M.print();
#endif
		Poly M0 = M.PXpPmX();
		M0 /= 2;
		M0 %= N;
		M0 = M0.sqSubst();
#ifdef VERBOSE
		cout << "M0: ";

		M0.print();
#endif
		Poly M1 = M.PXmPmX();
		M1 = (M1 >> 1);
		M1 /= 2;
		M1 %= N;
		M1 = M1.sqSubst();
#ifdef VERBOSE
		cout << "M1: ";
		M1.print();
#endif
		Poly M0sq = M0 * M0;
		Poly M1sq = M1 * M1;
		Poly Xpoly = Poly(1);
		Xpoly.set_coeff(1, 1);
		Poly V = M - M0sq;
		// V %= (1 << (N - newN));
		V = V + (Xpoly * M1sq);
		V /= (1 << newN);
		V %= (N - newN);
#ifdef VERBOSE
		cout << "V: ";
		V.print();
#endif
		Poly delta = teichmuller_modulus_increment(M0, M1, V, N - newN);
#ifdef VERBOSE
		cout << "d: ";
		delta.print();
#endif
		delta *= (1 << newN);

		M = M + delta;
		M %= N;
#ifdef VERBOSE
		cout << "M: ";
		M.print();
#endif
	}
	return M;
}

Poly Adicops::teichmuller_modulus_increment(const Poly& M0
		, const Poly& M1, const Poly& V, int N) {
	Poly delta;
	if (N == 1) {
		Poly nu(1);
		delta = nu - V;
		delta %= 1;
	} else {
		int newN = (N % 2 == 0 ? N / 2 : N / 2 + 1);
		delta = teichmuller_modulus_increment(M0, M1, V, newN);

		Poly delta0 = delta.PXpPmX();
		delta0 /= 2;
		delta0 %= N;
		delta0 = delta0.sqSubst();

		Poly delta1 = delta.PXmPmX();
		delta1 = (delta1 >> 1);
		delta1 /= 2;
		delta1 %= N;
		delta1 = delta1.sqSubst();

		Poly XPoly = Poly(1);
		XPoly.set_coeff(1, 1);
		Poly tmp = (delta0*M0) - (XPoly*M1*delta1);
		tmp *= 2;
		Poly Vnew = delta + V - tmp;
		Vnew /= (1 << newN);
		Vnew %= (N - newN);

		Poly bigDelta = teichmuller_modulus_increment(M0, M1, Vnew, N - newN);
		bigDelta *= (1 << newN);

		delta = bigDelta + delta;
		delta %= N;
	}

	return delta;
}

Poly Adicops::poly_invert(Poly f, int N, int prec) {
	if (N == 1) {
		return Poly::one();
	} else {
		int newN = (N % 2 == 0 ? N / 2 : N / 2 + 1);
		Poly c = poly_invert(f, newN, prec);
		Poly one = Poly::one();
		//Poly fc = c*(one - (f * c));
		c = c + (c*(one - (f * c)));
		//cout << "N: " << N << endl;
		//c.print();
		c = c.modXpowm(N);
		c %= prec;
		//c.print();
		return c;
	}
}

Poly Adicops::poly_remainder(Poly a, int prec) {
	return poly_remainder(a, mod, prec);
}

Poly Adicops::poly_remainder(Poly a, Poly b, int prec) {
	if (a.degree < b.degree) {
		a %= prec;
		return a;
	} else {
		int n = a.degree - b.degree + 1;
		//Poly powpoly = Poly(b.degree);
		//powpoly.set_coeff(b.degree, 1);
		Poly c = poly_invert(b.reverse(), n, prec);
		//c.print();
		Poly qq = a.reverse() * c;
		qq %= prec;
		qq = qq.modXpowm(n);
		Poly q = qq.reverse();
		Poly r = (a - (b*q));
		//cout << "r before, "; r.print();
		r %= prec;
		//cout << "r after, "; r.print();
		r = r.modXpowm(b.degree);
		return r;
	}
}

Poly Adicops::poly_division_rem(Poly a, Poly b, int prec) {
	if (a.degree < b.degree) {
		return Poly::zero();
	} else {
		int n = a.degree - b.degree + 1;
		//Poly powpoly = Poly(b.degree);
		//powpoly.set_coeff(b.degree, 1);
		Poly c = poly_invert(b.reverse(), n, prec);
		//c.print();
		Poly qq = a.reverse() * c;
		qq %= prec;
		qq = qq.modXpowm(n);
		Poly q = qq.reverse();
		return q;
	}
}

Poly Adicops::poly_division(Poly num, Poly denom, int prec) {
	//cout << "inv.. "; denom.print();
	Poly invDenom = get_inverse(denom, prec);
	Poly div = num * invDenom;
	div = poly_remainder(div, prec);
	return div;
}

Poly Adicops::get_inverse(Poly a, int prec) {
	//cout << "inv: "; a.print();
	if (prec == 1) {
		GFE gfe = GFE(a.to_gfe_el(), mod.to_gfe_el());
		//gfe.print();
		GFE invGfe = !gfe;
		//invGfe.print();

		//(gfe * invGfe).print();
		//mpz_class asd = invGfe.get_element();
		//cout << asd << endl;
		Poly res = Poly(invGfe.get_element());
		//cout << "N = " << prec << ":"; res.print();
		return res;
	} else {
		Poly z = get_inverse(a,
				(prec % 2 == 0 ? prec / 2 : prec / 2 + 1));
		Poly one = Poly::one();
		z = z + z*(one - (a*z));
		z = poly_remainder(z, prec);
		//z = z.modXpowm(prec);
		//z %= (1 << prec);
		//cout << "N = " << prec << ":"; z.print();
		return z;
	}

}

ModPoly Adicops::get_inverse(ModPoly a, int prec) {
	if (prec == 1) {
		ModPoly::set_precision(prec);
		GFE gfe = GFE(a.to_gfe_el(), mod.to_gfe_el());
		GFE invGfe = !gfe;
		ModPoly res = ModPoly(invGfe.get_element());
		return res;
	} else {
		ModPoly z = get_inverse(a,
				(prec % 2 == 0 ? prec / 2 : prec / 2 + 1));
		ModPoly::set_precision(prec);
		ModPoly one = ModPoly::one();
		cout << z.poly << endl;
		z = z + z * (one - (a * z));
		cout << z.poly << endl;
		return z;
	}
}

Poly Adicops::get_invsqrt(Poly a, Poly approx, int prec) {
	if (prec <= 2) {
		return approx;
	} else {
		int newN = ((prec + 1) % 2 == 0 ? (prec + 1) / 2 : (prec + 1) / 2 + 1);
		Poly z = get_invsqrt(a, approx, newN);

		Poly one = Poly::one();

		Poly amazsq = one - (a*z*z);
		amazsq = poly_remainder(amazsq, prec + 1);
		amazsq = z*amazsq;
		amazsq /= 2;
		z = z + amazsq;
		z = poly_remainder(z, prec);
		return z;
	}
}

ModPoly Adicops::get_invsqrt(ModPoly a, ModPoly approx, int prec) {
	if (prec <= 2) {
		return approx;
	} else {
		int newN = ((prec + 1) % 2 == 0 ? (prec + 1) / 2 : (prec + 1) / 2 + 1);
		ModPoly z = get_invsqrt(a, approx, newN);

		ModPoly::set_precision(prec + 1);
		ModPoly one = ModPoly::one();
		ModPoly amazsq = one - (a * z * z);
		//amazsq = poly_remainder(amazsq, prec + 1);
		amazsq = z * amazsq;
		amazsq /= 2;
		ModPoly::set_precision(prec);
		z = z + amazsq;
		//z = poly_remainder(z, prec);
		return z;
	}
}

Poly Adicops::get_sqrt(Poly a, int prec) {

	// binary inverse sqrt...
	GFE z = GFE(a.to_gfe_el(), mod.to_gfe_el());
	// TODO: ...
	GFE sqrtx = GFE::get_sqrtx(7, 1, mod.to_gfe_el());
	z = z.get_sqrt(sqrtx);
	z = !z;

	// Computing b=(1/a + z^2) / 4
	Poly inva = get_inverse(a, prec);
	Poly polyz = Poly(z.get_element());
	Poly b = inva - (polyz * polyz);
	b = poly_remainder(b, prec);
	b /= 4;

	// solving equation Delta^2+z*Delta=b
	GFE bgfe = GFE(b.to_gfe_el(), mod.to_gfe_el());
	GFE bigdelta = GFE::solve_quad_eq(z, bgfe);

	Poly polybigdelta = Poly(bigdelta.get_element());
	polybigdelta *= 2;

	// approx root to prec 2 is z+2*Delta
	polyz = polyz + polybigdelta;

	// comp. inverse square root with initial approximation polyz
	Poly invsqrt = get_invsqrt(a, polyz, prec);

	// revover sqrt: 1/a^{-1} * a = sqrt(a)
	invsqrt = invsqrt * a;
	return poly_remainder(invsqrt, prec);
}

ModPoly Adicops::get_sqrt(ModPoly a, int prec) {

	// binary inverse sqrt...
	GFE z = GFE(a.to_gfe_el(), mod.to_gfe_el());
	// TODO: ...
	GFE sqrtx = GFE::get_sqrtx(7, 1, mod.to_gfe_el());
	z = z.get_sqrt(sqrtx);
	z = !z;

	// Computing b=(1/a + z^2) / 4
	ModPoly inva = get_inverse(a, prec);
	ModPoly polyz = ModPoly(z.get_element());
	ModPoly::set_precision(prec);
	ModPoly b = inva - (polyz * polyz);
	//b = poly_remainder(b, prec);
	b /= 4;

	// solving equation Delta^2+z*Delta=b
	GFE bgfe = GFE(b.to_gfe_el(), mod.to_gfe_el());
	GFE bigdelta = GFE::solve_quad_eq(z, bgfe);

	ModPoly polybigdelta = ModPoly(bigdelta.get_element());
	polybigdelta *= 2;

	// approx root to prec 2 is z+2*Delta
	polyz = polyz + polybigdelta;

	// comp. inverse square root with initial approximation polyz
	ModPoly invsqrt = get_invsqrt(a, polyz, prec);
	// revover sqrt: 1/a^{-1} * a = sqrt(a)
	invsqrt = invsqrt * a;
	return invsqrt;
	//return poly_remainder(invsqrt, prec);
}

mpz_class Adicops::get_points_AGM_univariate(mpz_class _c, mpz_class _mod, int d) {
	int N = (d % 2 == 0 ? d / 2 + 3 : d / 2 + 4);
	Poly one = Poly::one();

	Poly c = Poly(_c);
	Poly mod = Poly(_mod);

	c *= 8;
	Poly b = one + c;
	Poly eps = poly_remainder(b, 4);
	cout << 4 << " eps: "; eps.print();

	for (int i = 5; i <= N; ++i) {
		Poly sqrtEps = get_sqrt(eps, i);
		assert(testsqrt(sqrtEps, eps, mod, i));
		//sqrtEps *= 2;
		//sqrtEps = poly_remainder(sqrtEps, mod, i);
		Poly t = one + eps;
		t /= 2;
		eps = poly_division(t, sqrtEps, i);
		//sqrtEps = poly_remainder(sqrtEps, mod, i);
		//sqrtEps = get_inverse(sqrtEps, mod, i);
		//eps = one + eps;
		//eps = eps * sqrtEps;
		//eps = poly_remainder(eps, mod, i);
		cout << i << " eps: "; eps.print();
	}

	Poly num = eps;
	num *= 2;
	Poly t = poly_division(num, one + eps, N - 1);
	return t.coeffs[0];
}

mpz_class Adicops::get_points_AGM_bivariate(mpz_class _c, int d) {

	int N = (d % 2 == 0 ? d / 2 + 3 : d / 2 + 4);
	Poly a = Poly::one();

	Poly c = Poly(_c);
	//Poly mod = Poly(_mod);

	c *= 8;
	Poly b = a + c;

	c = poly_remainder(b, 4);

	//cout << "4 a: "; a.print();
	//cout << "4 b: "; b.print();

	for (int i = 5; i <= N; ++i) {
		Poly olda = a;
		Poly oldb = b;
		a = a + b;
		a /= 2;
		Poly ab = olda * oldb;
		ab = poly_remainder(ab, i);
		b = get_sqrt(ab, i);
		//ab.print();
		assert(testsqrt(b, ab, mod, i));
		//Poly re = b * b;
		//re = poly_remainder(re, mod, i);
		cout << "N=" << i << " a: "; a.print();
		cout << "N=" << i << " b: "; b.print();
		//cout << i << " ab: "; ab.print();
		//cout << i << " re: "; re.print();
	}

	cout << "---" << endl;
	Poly a0 = a;

	for (int i = 0; i <= d - 1; ++i) {
		Poly olda = a;
		Poly oldb = b;
		a = olda + oldb;
		a /= 2;
		//a %= (1 << )
		Poly ab = olda * oldb;
		ab = poly_remainder(ab, N);
		b = get_sqrt(ab, N);
		Poly re = b * b;
		re = poly_remainder(re, N);
		assert(testsqrt(b, ab, mod, N));
		cout << "N=" << i << " a: "; a.print();
		cout << "N=" << i << " b: "; b.print();
		//cout << i << " ab: "; ab.print();
		//cout << i << " re: "; re.print();
	}

	a0 %= (N - 1);
	a %= (N - 1);


/*
	a0.set_coeff(0, 5);a0.set_coeff(1, 36);a0.set_coeff(2, 16);
	a0.set_coeff(3, 8);a0.set_coeff(4, 32);a0.set_coeff(5, 0);
	a0.set_coeff(6, 16);

	a.set_coeff(0, 57);a.set_coeff(1, 52);a.set_coeff(2, 48);
	a.set_coeff(3, 40);a.set_coeff(4, 32);a.set_coeff(5, 32);
	a.set_coeff(6, 48);
*/

	Poly t = poly_division(a0, a, N - 1);
	//Poly rem = poly_remainder(a0, a, N - 2);

	cout << "a0: "; a0.print();
	cout << "a: "; a.print();

	cout << "t: "; t.print();
	//cout << "rem: "; rem.print();

	//Poly inva = get_inverse(a, mod, N - 1);
	//Poly t = a0 * inva;
	//t = poly_remainder(t, mod, N - 1);
	//Poly t = poly_division(a0, a, N - 2);

	mpz_class c1 = 1;
	mpz_class mt = t.coeffs[0];
	//mpz_class mt = t.to_gfe_el();
	if ((mt * mt) > (c1 << (d + 2))) {
		mt = mt - (c1 << (N - 1));
	}

	return (c1 << d) + 1 - mt;
	//mod.print();
}

ZZ Adicops::get_points_AGM_bivariate_v2(mpz_class _c, int d) {
	int N = (d % 2 == 0 ? d / 2 + 3 : d / 2 + 4);
	ModPoly a = ModPoly::one();

	ModPoly c = ModPoly(_c);
	//cout << "modulus: " << ZZ_pE::modulus() << endl;
	//Poly mod = Poly(_mod);

	ModPoly::set_precision(4);
	c *= 8;
	ModPoly b = a + c;

	//c = poly_remainder(b, 4);

	//cout << "4 a: "; a.print();
	//cout << "4 b: "; b.print();

	for (int i = 5; i <= N; ++i) {
		ModPoly olda = a;
		ModPoly oldb = b;
		ModPoly::set_precision(i + 1);
		a = a + b;
		ModPoly::set_precision(i);
		cout << "a: " << a.poly << endl << "b: " << b.poly << endl;
		a /= 2;
		cout << "a: " << a.poly << endl << "b: " << b.poly << endl;
		ModPoly ab = olda * oldb;
		//ab = poly_remainder(ab, i);
		b = get_sqrt(ab, i);
		//ab.print();
		//assert(testsqrt(b, ab, mod, i));
		//Poly re = b * b;
		//re = poly_remainder(re, mod, i);
		cout << "N=" << i << " a: " << a.poly << endl;
		//a.print();
		cout << "N=" << i << " b: " << b.poly << endl;
		//b.print();
		//cout << i << " ab: "; ab.print();
		//cout << i << " re: "; re.print();
	}

	cout << "---" << endl;
	ModPoly a0 = a;

	for (int i = 0; i <= d - 1; ++i) {
		ModPoly olda = a;
		ModPoly oldb = b;
		ModPoly::set_precision(N + 1);
		a = olda + oldb;
		a /= 2;
		ModPoly::set_precision(N);
		//a %= (1 << )
		ModPoly ab = olda * oldb;
		//ab = poly_remainder(ab, N);
		b = get_sqrt(ab, N);
		ModPoly re = b * b;
		//re = poly_remainder(re, N);
		//assert(testsqrt(b, ab, mod, N));
		cout << "N=" << i << " a: " << a.poly << endl;
		//a.print();
		cout << "N=" << i << " b: " << b.poly << endl;
		//b.print();
		//cout << i << " ab: "; ab.print();
		//cout << i << " re: "; re.print();
	}

	//a0 %= (N - 1);
	//a %= (N - 1);

	ModPoly::set_precision(N - 1);
	//Poly t = poly_division(a0, a, N - 1);
	ModPoly t = a0 / a;
	//Poly rem = poly_remainder(a0, a, N - 2);

	cout << "a0: " << a0.poly << endl;
	//a0.print();
	cout << "a: " << a.poly << endl;
	//a.print();

	cout << "t: " << t.poly << endl;
	//t.print();
	//cout << "rem: "; rem.print();

	//Poly inva = get_inverse(a, mod, N - 1);
	//Poly t = a0 * inva;
	//t = poly_remainder(t, mod, N - 1);
	//Poly t = poly_division(a0, a, N - 2);

	//mpz_class c1 = 1;
	ZZ mt = t.poly.rep[0].cardinality();
	//mpz_class mt = t.to_gfe_el();
	if ((mt * mt) > (1 << (d + 2))) {
		mt = mt - (1 << (N - 1));
	}

	return (1 << d) + 1 - mt;
}

bool Adicops::testsqrt(Poly sqrt, Poly in, Poly mod, int prec) {
	Poly t = sqrt * sqrt;
	t = poly_remainder(t, prec);

	if (!(t == in)) {
		cout << "t1: ";
		t.print();
		cout << "t2: ";
		in.print();
	}
	return t == in;
}
