/*
 * ellipticcurve.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: bhess
 */

#include "ellipticcurve.h"

Ellipticcurve::Ellipticcurve() {}

Ellipticcurve::~Ellipticcurve() {}

Jacobian Ellipticcurve::addition(Jacobian P, Coordinate Q) {
    //I treat a**2 as simply a*a since I can't see an optimization
    //for it in gmp
    mpz_class X3,Y3,Z3;
    if (ECC_a != -3){
        //pg. 89 - dealing with the general case where a != -3
        mpz_class A,B,C,D,E,F,G,H,I;
        A = P.Z * P.Z;
        B = P.Z * A;
        C = Q.X * A;
        D = Q.Y * B;
        E = C - P.X;
        F = D - P.Y;
        G = E * E;
        H = G * E;
        I = P.X * G;
        X3 = F*F - (H+2*I);
        Y3 = F*(I-X3)-P.Y*H;
        Z3 = P.Z * E;
    }
    else{
        //pg. 91 - dealing with the case a == 3
        mpz_class T1,T2,T3,T4;
        T1 = P.Z * P.Z;
        T2 = T1 * P.Z;
        T1 = T1 * Q.X;
        T2 = T2 * Q.Y;
        T1 = T1 - P.X;
        T2 = T2 - P.Y;
        if (T1 == 0)
            if (T2 == 0)
                return doubling(Q);
            else
                return Jacobian(1,1,0); //TODO: write a class method to return the identity
        Z3 = P.Z * T1;
        T3 = T1 * T1;
        T4 = T3 * T1;
        T3 = T3 * P.X;
        T1 = 2 * T3;
        X3 = T2 * T2;
        X3 = X3 - T1;
        X3 = X3 - T4;
        T3 = T3 - X3;
        T3 = T3 * T2;
        T4 = T4 * P.Y;
        Y3 = T3 - T4;
    }
    X3 = X3 % mod;
    Y3 = Y3 % mod;
    Z3 = Z3 % mod;
    return Jacobian(X3,Y3,Z3);
}

Jacobian Ellipticcurve::doubling(Jacobian P) {
    //I treat a**2 as simply a*a since I can't see an optimization
    //for it in gmp
    mpz_class X3,Y3,Z3;
    if (ECC_a != -3){
        //pg. 88 - dealing with the general case where a != -3
        mpz_class A,B,C,D;
        A = P.Y * P.Y;
        B = 4*P.X * A;
        C = 8 * A * A;
        D = 3*P.X*P.X + ECC_a*P.Z*P.Z*P.Z*P.Z;
        X3 = D*D-2*B;
        Y3 = D*(B-X3) - C;
        Z3 = 2*P.Y*P.Z;
    }
    else{
        //pg. 91 - dealing with the case a==-3
        mpz_class T1, T2, T3;
        T1 = P.Z*P.Z;
        T2 = P.X - T1;
        T1 = P.X + T1;
        T2 = T2 * T1;
        T2 = 3 * T2;
        Y3 = 2 * P.Y;
        Z3 = Y3 * P.Z;
        Y3 = Y3 * Y3;
        T3 = Y3 * P.X;
        Y3 = Y3 * Y3;
        Y3 = Y3 / 2;
        X3 = T2 * T2;
        T1 = 2 * T3;
        X3 = X3 - T1;
        T1 = T3 - X3;
        T1 = T1 * T2;
        Y3 = T1 - Y3;
    }
    X3 = X3 % mod;
    Y3 = Y3 % mod;
    Z3 = Z3 % mod;
    return Jacobian(X3,Y3,Z3);
}

Coordinate Ellipticcurve::pointMultiplication(Coordinate P, mpz_class k) {

}

Jacobian Ellipticcurve::repeatedDoubling(Jacobian P, int m) {
	// For the case a==-3
	// TODO: test...

	if (P.isInfinite()) {
		return P;
	} else {
		mpz_class Y, W, A, B, X, Z, tmp_W;
		X = P.X; Y = P.Y; Z = P.Z;
		mpz_class const_4 = 4;

		Y *= 2;
		mpz_powm(W.get_mpz_t(), P.Z.get_mpz_t(), const_4.get_mpz_t(), mod.get_mpz_t());
		while (m > 0) {
			A = 3 * (X*X - W);
			B = X * Y * Y;
			X = A*A - 2*B;
			Z = Z * Y;
			if (--m > 0) {
				tmp_W = W;
				mpz_powm(W.get_mpz_t(), Y.get_mpz_t(), const_4.get_mpz_t(), mod.get_mpz_t());

				// avoid computing y^4 two times -> W is now Y^4
				Y = 2*A * (B - X) - W;

				W = tmp_W * W;
			} else {
				// Y <-mod Y^4
				mpz_powm(Y.get_mpz_t(), Y.get_mpz_t(), const_4.get_mpz_t(), mod.get_mpz_t());
				Y = 2*A * (B - X) - Y;
			}
		}

		X = X % mod;
		Y = (Y / 2) % mod;
		Z = Z % mod;

		return Jacobian(X, Y, Z);
	}
}

std::vector<int> Ellipticcurve::getNAF(mpz_class k) {
    //implementation folows pg. 98
    //TODO: test
    std::vector<int> naf;
    int i = 0;
    while (k >= 1){
        if (k % 2 == 1){
            mpz_class temp = (k % 4);
            naf[i] = 2 - temp.get_ui();
            k -= naf[i];
        }
        else{
            naf[i] = 0;
        }
        k /= 2;
        i++;
    }
    return naf;
}
