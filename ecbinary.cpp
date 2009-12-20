/*
 * ECBinary.cpp
 *
 *  Created on: Nov 22, 2009
 *      Author: bhess
 */

#include "ecbinary.h"

ECBinary::ECBinary() {
	// TODO Auto-generated constructor stub

}

Coordinate ECBinary::addition(Coordinate P, Coordinate Q) {
	return toCoordinate(addition(LD(P), Q));
}

Coordinate ECBinary::subtraction(Coordinate P, Coordinate Q) {
	return toCoordinate(subtraction(LD(P), Q));
}

Coordinate ECBinary::doubling(Coordinate P) {
	return toCoordinate(doubling(LD(P)));
}

Coordinate ECBinary::repeatedDoubling(Coordinate P, int m) {
	// not implemented...
}

Coordinate ECBinary::pointMultiplication(Coordinate P, mpz_class k) {

	// implementation according to p.99
	// actually same implementation as with prime curves,
	//  but using LD coordinates

	std::vector<int> naf = getNAF(k);
	LD Q = LD::infinity();
	for (int i = naf.size() - 1; i >= 0; --i) {
		Q = doubling(Q);
		if (naf[i] == 1) {
			Q = addition(Q, P);
		} else if (naf[i] == -1) {
			Q = subtraction(Q, P);
		}
	}
	return Coordinate(Q, mod);

}

/**
 * Optimization when a\in {0,1}
 * according p.94
 *
 * TODO: test
 */
LD ECBinary::addition(LD P, Coordinate Q) {

	// identity property
	if (Q.isInfinite()) {
		return P;
	}
	if (P.isInfinite()) {
		return LD(Q);
	}
	mpz_class x3, y3, z3;
	mpz_class const_2 = 2;
	if (ECC_a == 0 || ECC_a == 1) {
		mpz_class t1, t2, t3, y3;
		t1 = P.Z * Q.X;
		mpz_mod(t1.get_mpz_t(), t1.get_mpz_t(), mod.get_mpz_t());
		mpz_powm(t2.get_mpz_t(), P.Z.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		x3 = P.X + t1;
		mpz_mod(x3.get_mpz_t(), x3.get_mpz_t(), mod.get_mpz_t());
		t1 = P.Z * x3;
		mpz_mod(t1.get_mpz_t(), t1.get_mpz_t(), mod.get_mpz_t());
		t3 = t2 * Q.Y;
		mpz_mod(t3.get_mpz_t(), t3.get_mpz_t(), mod.get_mpz_t());
		y3 = P.Y + t3;
		mpz_mod(y3.get_mpz_t(), y3.get_mpz_t(), mod.get_mpz_t());
		if (x3 == 0) {
			return (y3 == 0 ? doubling(LD(Q)) : LD::infinity());
		}
		mpz_powm(z3.get_mpz_t(), t1.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		t3 = t1 * y3;
		mpz_mod(t3.get_mpz_t(), t3.get_mpz_t(), mod.get_mpz_t());
		if (ECC_a == 1) {
			t1 = t1 + t2;
			mpz_mod(t1.get_mpz_t(), t1.get_mpz_t(), mod.get_mpz_t());
		}
		mpz_powm(t2.get_mpz_t(), x3.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		x3 = t2 * t1;
		mpz_mod(x3.get_mpz_t(), x3.get_mpz_t(), mod.get_mpz_t());
		mpz_powm(t2.get_mpz_t(), y3.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		x3 = x3 + t2;
		mpz_mod(x3.get_mpz_t(), x3.get_mpz_t(), mod.get_mpz_t());
		x3 = x3 + t3;
		mpz_mod(x3.get_mpz_t(), x3.get_mpz_t(), mod.get_mpz_t());
		t2 = Q.X * z3;
		mpz_mod(t2.get_mpz_t(), t2.get_mpz_t(), mod.get_mpz_t());
		t2 = t2 + x3;
		mpz_mod(t2.get_mpz_t(), t2.get_mpz_t(), mod.get_mpz_t());
		mpz_powm(t1.get_mpz_t(), z3.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		t3 = t3 + z3;
		mpz_mod(t3.get_mpz_t(), t3.get_mpz_t(), mod.get_mpz_t());
		y3 = t3 * t2;
		mpz_mod(y3.get_mpz_t(), y3.get_mpz_t(), mod.get_mpz_t());
		t2 = Q.X + Q.Y;
		mpz_mod(t2.get_mpz_t(), t2.get_mpz_t(), mod.get_mpz_t());
		t3 = t1 * t2;
		mpz_mod(t3.get_mpz_t(), t3.get_mpz_t(), mod.get_mpz_t());
		y3 = y3 + t3;
		mpz_mod(y3.get_mpz_t(), y3.get_mpz_t(), mod.get_mpz_t());
	} else {
		mpz_class a, b, c, d, e, f, g;
		a = (Q.Y * P.Z * P.Z) + P.Y;
		mpz_mod(a.get_mpz_t(), a.get_mpz_t(), mod.get_mpz_t());
		b = Q.X * P.Z + P.X;
		mpz_mod(b.get_mpz_t(), b.get_mpz_t(), mod.get_mpz_t());
		c = P.Z * b;
		mpz_mod(c.get_mpz_t(), c.get_mpz_t(), mod.get_mpz_t());
		d = b * b * (c + ECC_a * P.Z * P.Z);
		mpz_mod(d.get_mpz_t(), d.get_mpz_t(), mod.get_mpz_t());
		mpz_powm(z3.get_mpz_t(), c.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		e = a * c;
		mpz_mod(e.get_mpz_t(), e.get_mpz_t(), mod.get_mpz_t());
		x3 = a * a + d + e;
		mpz_mod(x3.get_mpz_t(), x3.get_mpz_t(), mod.get_mpz_t());
		f = x3 + Q.X * z3;
		mpz_mod(f.get_mpz_t(), f.get_mpz_t(), mod.get_mpz_t());
		g = (Q.X + Q.Y) * z3 * z3;
		mpz_mod(g.get_mpz_t(), g.get_mpz_t(), mod.get_mpz_t());
		y3 = (e + z3) * f + g;
		mpz_mod(y3.get_mpz_t(), y3.get_mpz_t(), mod.get_mpz_t());
	}
	return LD(x3, y3, z3);
}

LD ECBinary::subtraction(LD P, Coordinate Q) {
	return addition(P, getNegative(Q));
}

/*
 * Optimized for the case when a=1 or 2
 * according to p.94
 *
 * TODO: test...
 */
LD ECBinary::doubling(LD P) {
	if (P.isInfinite()) {
		return P;
	}

	mpz_class x3, y3, z3;
	mpz_class const_2 = 2;
	if (ECC_a == 0 || ECC_a == 1) {
		mpz_class t1, t2;
		mpz_powm(t1.get_mpz_t(), P.Z.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		mpz_powm(t2.get_mpz_t(), P.X.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		z3 = t1 * t2;

		mpz_powm(x3.get_mpz_t(), t2.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		mpz_powm(t1.get_mpz_t(), t1.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		t2 = t1 * ECC_b;
		mpz_mod(t2.get_mpz_t(), t2.get_mpz_t(), mod.get_mpz_t());
		x3 = x3 + t2;
		mpz_mod(x3.get_mpz_t(), x3.get_mpz_t(), mod.get_mpz_t());
		mpz_powm(t1.get_mpz_t(), P.Y.get_mpz_t(), const_2.get_mpz_t(),
				mod.get_mpz_t());
		if (ECC_a == 1) {
			t1 = t1 + z3;
			mpz_mod(t1.get_mpz_t(), t1.get_mpz_t(), mod.get_mpz_t());
		}
		t1 = t1 + t2;
		mpz_mod(t1.get_mpz_t(), t1.get_mpz_t(), mod.get_mpz_t());
		y3 = x3 * t1;
		mpz_mod(y3.get_mpz_t(), y3.get_mpz_t(), mod.get_mpz_t());
		t1 = t2 * z3;
		mpz_mod(t1.get_mpz_t(), t1.get_mpz_t(), mod.get_mpz_t());
		y3 = y3 + t1;
		mpz_mod(y3.get_mpz_t(), y3.get_mpz_t(), mod.get_mpz_t());
	} else {
		// TODO: intermediate results...
		z3 = P.X * P.X * P.Z * P.Z;
		mpz_mod(z3.get_mpz_t(), z3.get_mpz_t(), mod.get_mpz_t());
		x3 = P.X * P.X * P.X * P.X + ECC_b * P.Z * P.Z * P.Z * P.Z;
		mpz_mod(x3.get_mpz_t(), x3.get_mpz_t(), mod.get_mpz_t());
		y3 = ECC_b * P.Z * P.Z * P.Z * P.Z * z3 + x3 * (ECC_a * z3 + P.Y * P.Y
				+ ECC_b * P.Z * P.Z * P.Z * P.Z);
		mpz_mod(y3.get_mpz_t(), y3.get_mpz_t(), mod.get_mpz_t());
	}
	return LD(x3, y3, z3);
}

Coordinate ECBinary::toCoordinate(const LD& ld) {

	GFE ldX = GFE(ld.X, mod), ldY = GFE(ld.Y, mod), ldZ = GFE(ld.Z, mod);
	GFE resY = ldZ * ldZ;
	GFE resX = !ldZ;
	resY = !resY;
	resX = ldX * resX;
	resY = ldY * resY;

	return Coordinate(resX.element, resY.element);

}

vector<GFE> ECBinary::get_jinvariants() {
	GFE a = GFE(ECC_b, mod);
	a.print();
	// j(E)=j0=1/b
	GFE j0 = !a;
	// We create cyclic curves, therefore: jd=j0
	GFE jd = j0;
	int d = a.mod_deg();
	cout << d << endl;
	std::vector<GFE> res(d + 1);
	//res[0] = j0;
	res[d] = jd;
	for (int i = d - 1; i >= 0; --i) {
		res[i] = res[i + 1] * res[i + 1];
	}

	return res;
}


std::vector<GFE> ECBinary::update_js(int n, std::vector<GFE> Jinvs) {
	int d = Jinvs.size();
	vector<GFE> D(d - 1);
	vector<GFE> P(d);
	vector<GFE> J(d);
	for (int i = 0; i < d - 2; ++i) {
		GFE t = !phi_2_x(Jinvs[i], Jinvs[i+1]);

	}


	return J;
}

/*
 * Lifts j-invariants given by 'jinvs' to precision n
 * jinvs: j_0,...,j_d
 * returns J_0,...,J_{d-1}
 */
std::vector<GFE> ECBinary::lift_jinvariants(int n, std::vector<GFE> jinvs) {
	if (n == 1) {
		return jinvs;
	}
	int nn = (n % 2 == 0 ? n / 2 : n / 2 + 1);
	vector<GFE> res = lift_jinvariants(nn, jinvs);
	return update_js(n, res);
	// phi_2(X,Y)=X^3+Y^3-X^2Y^2+c1(XY^2+X^2Y)-c2(X^2+Y^2)+c3XY+c4(X+Y)-c5
	// c1=1488 c2=162000 c3=40773375 c4=87448000000 c5=157464000000000
}

mpz_class ECBinary::satohfgh_point_counting() {
	vector<GFE> jinvs = get_jinvariants();
	int d = jinvs[0].mod_deg();
	int n = (d % 2 == 0 ? d / 2 + 1 : d / 2 + 2);

	cout << "n: " << n << endl;
	return n;
}

inline GFE ECBinary::phi_2(GFE x, GFE y) {
	mpz_class c4t;
	c4t.set_str("87448000000", 2);
	mpz_class c5t;
	c5t.set_str("157464000000000", 2);
	GFE c1 = GFE(1488, x.mod);
	GFE c2 = GFE(162000, x.mod);
	GFE c3 = GFE(40773375, x.mod);
	GFE c4 = GFE(c4t, x.mod);
	GFE c5 = GFE(c5t, x.mod);
	return x*x*x+y*y*y-x*x*y*y+c1*(x*y*y+x*x*y)-c2*(x*x+y*y)+c3*x*y+c4*(x+y)-c5;
}

inline GFE ECBinary::phi_2_x(GFE x, GFE y) {
	GFE co2 = GFE(2, x.mod);
	GFE co3 = GFE(3, x.mod);
	GFE c1 = GFE(1488, x.mod);
	GFE c2 = GFE(162000, x.mod);
	GFE c3 = GFE(40773375, x.mod);
	return co3*x*x - co2*x*(y*y - c1*y + c2) + c1*y*y + c3*y;
}
