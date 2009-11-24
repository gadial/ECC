/*
 * ECPrime.cpp
 *
 *  Created on: Nov 22, 2009
 *      Author: bhess
 */

#include "ecprime.h"

ECPrime::ECPrime() {
	// TODO Auto-generated constructor stub

}

ECPrime::~ECPrime() {
	// TODO Auto-generated destructor stub
}

ECPrime ECPrime::randomCurve(int number_of_bits, RandomNumberGenerator gen){
    mpz_class p = gen.generate_prime(number_of_bits);
    mpz_class a = gen.rand(p);
    mpz_class b = gen.rand(p);

    //TODO: check that the curve is smooth

    return ECPrime(p,a,b);
}

Coordinate ECPrime::addition(Coordinate P, Coordinate Q) {
	return Coordinate(addition(Jacobian(P), Q), mod);
}

Coordinate ECPrime::subtraction(Coordinate P, Coordinate Q) {
	return Coordinate(subtraction(Jacobian(P), Q), mod);
}

Coordinate ECPrime::doubling(Coordinate P) {
	return Coordinate(doubling(Jacobian(P)), mod);
}

Coordinate ECPrime::repeatedDoubling(Coordinate P, int m) {
	return Coordinate(repeatedDoubling(Jacobian(P), m), mod);
}

Jacobian ECPrime::addition(Jacobian P, Coordinate Q) {
    //We assume that P != Q,-Q
    //If P==Q use doubling
    //if P==-Q, return 0
    //I treat a**2 as simply a*a since I can't see an optimization
    //for it in gmp
	if (P.isInfinite()) {
		return Jacobian(Q);
	}
	if (Q.isInfinite()) {
		return P;
	}
        Coordinate temp = Coordinate(P,mod);
        if (temp == Q)
            return doubling(P);
        if (temp == getNegative(Q))
            return Coordinate::infinity();
        
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
        X3 = F*F - (H+2*I); //cout << "F = " << F << ", H= " << H << ", I=" << I << endl;
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
                return doubling(Jacobian(Q));
            else
                return Coordinate::infinity();
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
    mpz_mod(X3.get_mpz_t(), X3.get_mpz_t(), mod.get_mpz_t());
    mpz_mod(Y3.get_mpz_t(), Y3.get_mpz_t(), mod.get_mpz_t());
    mpz_mod(Z3.get_mpz_t(), Z3.get_mpz_t(), mod.get_mpz_t());
    //X3 = X3 % mod;
    //Y3 = Y3 % mod;
    //Z3 = Z3 % mod;
    return Jacobian(X3,Y3,Z3);
}

Jacobian ECPrime::subtraction(Jacobian P, Coordinate Q) {
	return addition(P, getNegative(Q));
}

Jacobian ECPrime::doubling(Jacobian P) {
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

Coordinate ECPrime::pointMultiplication(Coordinate P, mpz_class k) {

	// implementation according to p.99

	std::vector<int> naf = getNAF(k);
	Jacobian Q = Jacobian::infinity();
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

Jacobian ECPrime::repeatedDoubling(Jacobian P, int m) {
	// TODO: test...

	if (P.isInfinite()) {
		return P;
        } else if (ECC_a != -3){ //naive, for the case a != -3
            Jacobian temp = P;
            for (int i=1; i<=m; i++)
                temp = doubling(temp);
            return temp;
	} else { // For the case a==-3
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
