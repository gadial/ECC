/*
 * ellipticcurve.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: bhess
 */

#include <iostream>

#include "primes.h"

#include "ellipticcurve.h"

using namespace std;
//#include <cassert>

Ellipticcurve::Ellipticcurve() {}

Ellipticcurve::~Ellipticcurve() {}

Ellipticcurve::Ellipticcurve(const char* _mod, int _mod_base,
			const char* _order, int _order_base,
			const char* _ecc_a, int _ecc_a_base,
			const char* _ecc_b, int _ecc_b_base,
			const char* _px, int _px_base,
			const char* _py, int _py_base) {

	mod.set_str(_mod, _mod_base);
	order.set_str(_order, _order_base);
	ECC_a.set_str(_ecc_a, _ecc_a_base);
	ECC_b.set_str(_ecc_b, _ecc_b_base);
	point = Coordinate(_px, _px_base, _py, _py_base);
}

std::vector<int> Ellipticcurve::getNAF(mpz_class k) {
    //implementation folows pg. 98
    std::vector<int> naf(mpz_sizeinbase(k.get_mpz_t(), 2) + 1);
    int i = 0;
    while (k >= 1){
        if (k % 2 == 1){
            mpz_class temp = (k % 4);
            naf[i] = 2 - (int)temp.get_ui();
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

Coordinate Ellipticcurve::getNegative(const Coordinate& P) {
	return Coordinate(P.X, mod - P.Y);
}

Coordinate Ellipticcurve::getPoint(mpz_class x, bool negative_value)
{
    //we solve the equation y^2 = x^3+ax+b
    mpz_class temp = (x*x*x + ECC_a*x + ECC_b) % mod;
    mpz_class y = modular_square_root(temp,mod);
    if (y == 0)
        return Coordinate::infinity();
    if (negative_value)
        return Coordinate(x,mod-y);
    return Coordinate(x,y);
}
