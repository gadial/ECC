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
	return Coordinate(addition(LD(P), Q), mod);
}

Coordinate ECBinary::subtraction(Coordinate P, Coordinate Q) {
	return Coordinate(subtraction(LD(P), Q), mod);
}

Coordinate ECBinary::doubling(Coordinate P) {
	return Coordinate(doubling(LD(P)), mod);
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
