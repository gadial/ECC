/*
 * ECPrime.cpp
 *
 *  Created on: Nov 22, 2009
 *      Author: bhess
 */

#include "ecprime.h"
#include "hcp.h"

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
bool ECPrime::check_order(mpz_class order_candidate){
    RandomNumberGenerator gen;
    for (int i=0; i<10; i++){
        zp_int x = gen.generate_modulu_p(mod);
        ZpCoordinate P = getPoint(x);
        
    }
}
static ECPrime::ECPrime randomCurveFromDiscriminant(int D, int number_of_bits, RandomNumberGenerator gen){
    #define SMOOTHNESS_ALLOWED 2
    #define MIN_SIZE_ALLOWED 10000
    mpz_class t,s;
    mpz_class u1,u2;
    mpz_class p;
    mpz_class u = 0;
    while (u == 0){
        p = gen.generate_prime_for_discriminant(number_of_bits,D, t,s); //4p=t^2+Ds^2
        u1 = p+1-t;
        u2 = p+1+t;

        if (is_near_prime(u1,SMOOTHNESS_ALLOWED,MIN_SIZE_ALLOWED))
            u = u1;
        if (is_near_prime(u2,SMOOTHNESS_ALLOWED,MIN_SIZE_ALLOWED))
            u = u2;
    }
    ModularPolynomial pol = ModularPolynomial::build_hcp_from_discriminant(D,p);
    zp_int j0 = pol.find_one_root();
    zp_int k = j0/(zp_int(1728,p)-j0);
    zp_int c = gen.generate_modulu_p(p);
    zp_int a = k*3*(c^2);
    zp_int b = k*(c^3);

    ECPrime candidate(p,a,b);

}

ZpCoordinate ECPrime::addition(ZpCoordinate P, ZpCoordinate Q) {
	return ZpCoordinate(addition(ZpJacobian(P), Q));
}

ZpCoordinate ECPrime::subtraction(ZpCoordinate P, ZpCoordinate Q) {
	return ZpCoordinate(subtraction(ZpJacobian(P), Q));
}

ZpCoordinate ECPrime::doubling(ZpCoordinate P) {
//    cout << "(doubling) P.p = " << P.p << endl;
//    cout << "(doubling) jac(P).p = " << ZpJacobian(P).Z.get_p() << endl;
//    cout << "(doubling) d(jac(P)).p = " << doubling(ZpJacobian(P)).Z.get_p() << endl;
//    cout << "(doubling) output.p = " << ZpCoordinate(doubling(ZpJacobian(P))).p << endl;
	return ZpCoordinate(doubling(ZpJacobian(P)));
}

ZpCoordinate ECPrime::repeatedDoubling(ZpCoordinate P, int m) {
	return ZpCoordinate(repeatedDoubling(ZpJacobian(P), m));
}

ZpJacobian ECPrime::addition(ZpJacobian P, ZpCoordinate Q) {
    //We assume that P != Q,-Q
    //If P==Q use doubling
    //if P==-Q, return 0
    //I treat a**2 as simply a*a since I can't see an optimization
    //for it in gmp
	if (P.isInfinite()) {
		return ZpJacobian(Q);
	}
	if (Q.isInfinite()) {
		return P;
	}
        ZpCoordinate temp = ZpCoordinate(P);
        if (temp == Q){
            return doubling(P);
        }
        if (temp == getNegative(Q)){
            return ZpJacobian::infinity(P.p);
        }
        
    zp_int X3,Y3,Z3;
    if (ECC_a != -3){
        //pg. 89 - dealing with the general case where a != -3
        zp_int A,B,C,D,E,F,G,H,I;
        A = P.Z * P.Z;
        B = P.Z * A;
        C = Q.X * A;
        D = Q.Y * B;
        E = C - P.X;
        F = D - P.Y;
        G = E * E;
        H = G * E;
        I = P.X * G;
        X3 = F*F - (H+I*2); //cout << "F = " << F << ", H= " << H << ", I=" << I << endl;
        Y3 = F*(I-X3)-P.Y*H;
        Z3 = P.Z * E;
    }
    else{
        //pg. 91 - dealing with the case a == 3
        zp_int T1,T2,T3,T4;
        T1 = P.Z * P.Z;
        T2 = T1 * P.Z;
        T1 = T1 * Q.X;
        T2 = T2 * Q.Y;
        T1 = T1 - P.X;
        T2 = T2 - P.Y;
        if (T1 == 0)
            if (T2 == 0)
                return doubling(ZpJacobian(Q));
            else
                return ZpCoordinate::infinity();
        Z3 = P.Z * T1;
        T3 = T1 * T1;
        T4 = T3 * T1;
        T3 = T3 * P.X;
        T1 = T3 * 2;
        X3 = T2 * T2;
        X3 = X3 - T1;
        X3 = X3 - T4;
        T3 = T3 - X3;
        T3 = T3 * T2;
        T4 = T4 * P.Y;
        Y3 = T3 - T4;
    }
    return ZpJacobian(X3,Y3,Z3);
}

ZpJacobian ECPrime::subtraction(ZpJacobian P, ZpCoordinate Q) {
	return addition(P, getNegative(Q));
}

ZpJacobian ECPrime::doubling(ZpJacobian P) {
    //I treat a**2 as simply a*a since I can't see an optimization
    //for it in gmp
    if (P.isInfinite())
        return P;
//    cout << "P = " << P << endl;
    zp_int X3,Y3,Z3;
    if (ECC_a != -3){
        //pg. 88 - dealing with the general case where a != -3
        zp_int A,B,C,D;
        A = P.Y * P.Y;
        B = P.X*4 * A;
        C = A*8 * A;
        D = P.X*P.X*3 + P.Z*P.Z*P.Z*P.Z*ECC_a;
        X3 = D*D-(B*2);
        Y3 = D*(B-X3) - C;
        Z3 = P.Y*P.Z*2;
    }
    else{
        //pg. 91 - dealing with the case a==-3
        zp_int T1, T2, T3;
        T1 = P.Z*P.Z;
        T2 = P.X - T1;
        T1 = P.X + T1;
        T2 = T2 * T1;
        T2 = T2 * 3;
        Y3 = P.Y * 2;
        Z3 = Y3 * P.Z;
        Y3 = Y3 * Y3;
        T3 = Y3 * P.X;
        Y3 = Y3 * Y3;
        Y3 = Y3 / 2;
        X3 = T2 * T2;
        T1 = T3 * 2;
        X3 = X3 - T1;
        T1 = T3 - X3;
        T1 = T1 * T2;
        Y3 = T1 - Y3;
    }
//    X3 = X3 % mod;
//    Y3 = Y3 % mod;
//    Z3 = Z3 % mod;
    return ZpJacobian(X3,Y3,Z3);
}

ZpCoordinate ECPrime::pointMultiplication(ZpCoordinate P, mpz_class k) {

	// implementation according to p.99

	std::vector<int> naf = getNAF(k);
	ZpJacobian Q = ZpJacobian::infinity(P.p);
	for (int i = naf.size() - 1; i >= 0; --i) {
		Q = doubling(Q);
		if (naf[i] == 1) {
			Q = addition(Q, P);
		} else if (naf[i] == -1) {
			Q = subtraction(Q, P);
		}
	}
	return ZpCoordinate(Q);
}

ZpJacobian ECPrime::repeatedDoubling(ZpJacobian P, int m) {
	// TODO: test...

	if (P.isInfinite()) {
		return P;
        } else if (ECC_a != -3){ //naive, for the case a != -3
            ZpJacobian temp = P;
            for (int i=1; i<=m; i++)
                temp = doubling(temp);
            return temp;
	} else { // For the case a==-3
		zp_int Y, W, A, B, X, Z, tmp_W;
		X = P.X; Y = P.Y; Z = P.Z;
		mpz_class const_4 = 4;

		Y *= 2;
                W = (P.Z^const_4);
		while (m > 0) {
			A = (X*X - W) * 3;
			B = X * Y * Y;
			X = A*A - B*2;
			Z = Z * Y;
			if (--m > 0) {
				tmp_W = W;
                                W = (Y^const_4);
				// avoid computing y^4 two times -> W is now Y^4
				Y = A * (B - X) * 2 - W;

				W = tmp_W * W;
			} else {
				// Y <-mod Y^4
                                Y = (Y^const_4);
				Y = A * (B - X) * 2 - Y;
			}
		}

		Y = (Y / 2);

		return ZpJacobian(X, Y, Z);
	}
}

ZpCoordinate ECPrime::getPoint(zp_int x, bool negative_value){
//    cout << "generating modulu x.p = " << x.get_p() << endl;
    zp_int y = modular_square_root(x*x*x + ECC_a*x + ECC_b);
//    cout << "generating modulu y.p = " << y.get_p() << endl;
    if (y == 0)
        return ZpCoordinate::infinity();
    if (negative_value)
        return ZpCoordinate(x,-y,mod);
    return ZpCoordinate(x,y,mod);
}

ZpCoordinate ECPrime::getNegative(const ZpCoordinate& P){
    return ZpCoordinate(P.X, -P.Y,P.p);
}