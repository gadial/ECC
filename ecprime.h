/*
 * ECPrime.h
 *
 *  Created on: Nov 22, 2009
 *      Author: bhess
 */

#ifndef ECPRIME_H_
#define ECPRIME_H_

#include "ellipticcurve.h"
#include "zp_int.h"
class ECPrime : public Ellipticcurve{
public:
	ECPrime();
	ECPrime(mpz_class _mod, mpz_class _ECC_a, mpz_class _ECC_b) :
            Ellipticcurve(_mod, _ECC_a, _ECC_b) {point = get_random_point();}
        ECPrime(const char* _mod, int _mod_base, const char* _order,
            int _order_base, const char* _ecc_a, int _ecc_a_base,
            const char* _ecc_b, int _ecc_b_base, const char* _px, int _px_base,
            const char* _py, int _py_base) :
                Ellipticcurve(_mod, _mod_base, _order, _order_base, _ecc_a, _ecc_a_base,_ecc_b, _ecc_b_base, _px, _px_base, _py, _py_base) {}
        
	static ECPrime randomCurve(int number_of_bits,
			RandomNumberGenerator gen);

        static ECPrime randomCurveFromDiscriminant(int D, int number_of_bits,
			RandomNumberGenerator gen);

        static ECPrime normalizedCurveFromDiscriminantAndPrime(int D, mpz_class p, int HCP_root_number = 0);
	/**
	 * Addition P+Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	ZpCoordinate addition(ZpCoordinate P, ZpCoordinate Q);

	/**
	 * Subtraction P-Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	ZpCoordinate subtraction(ZpCoordinate P, ZpCoordinate Q);

	/**
	 * Doubling of a point P -> 2P
	 * in jacobian coordinates
	 */
	ZpCoordinate doubling(ZpCoordinate P);

	/**
	 * Repeated doubling of a point P in jacobian coordinates
	 * (m times) -> 2^m P
	 */
	ZpCoordinate repeatedDoubling(ZpCoordinate P, int m);

	/**
	 * Point multiplication (k times)
	 * -> kP
	 */
	ZpCoordinate pointMultiplication(ZpCoordinate P, mpz_class k);

	virtual ~ECPrime();
        ZpCoordinate getPoint(zp_int x, bool negative_value = false);
        ZpCoordinate getPoint(mpz_class x, bool negative_value = false){return getPoint(zp_int(x,mod),negative_value);}
        ZpCoordinate getPointFromCompressedForm(string form);
        Coordinate addition(Coordinate P, Coordinate Q){return addition(ZpCoordinate(P,mod),ZpCoordinate(Q,mod));}
	Coordinate subtraction(Coordinate P, Coordinate Q){return subtraction(ZpCoordinate(P,mod),ZpCoordinate(Q,mod));}
	Coordinate doubling(Coordinate P){return doubling(ZpCoordinate(P,mod));}
	Coordinate repeatedDoubling(Coordinate P, int m){return repeatedDoubling(ZpCoordinate(P,mod),m);}
	Coordinate pointMultiplication(Coordinate P, mpz_class k){return pointMultiplication(ZpCoordinate(P,mod),k);}
	Coordinate getPoint_interface(mpz_class x, bool negative_value = false){return getPoint(x, negative_value);}
	Coordinate getPointCompressedForm(string from) {return getPointFromCompressedForm(from); }
        
        ZpCoordinate get_random_point();
        bool check_order(mpz_class order_candidate);
private:
	/**
	 * Addition P+Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	ZpJacobian addition(ZpJacobian P, ZpCoordinate Q);

	/**
	 * Subtraction P-Q of a jacobian coordinate
	 * P and an affine coordinate Q
	 */
	ZpJacobian subtraction(ZpJacobian P, ZpCoordinate Q);

	/**
	 * Doubling of a point P -> 2P
	 * in jacobian coordinates
	 */
	ZpJacobian doubling(ZpJacobian P);

	/**
	 * Repeated doubling of a point P in jacobian coordinates
	 * (m times) -> 2^m P
	 */
	ZpJacobian repeatedDoubling(ZpJacobian P, int m);

protected:
    ZpCoordinate getNegative(const ZpCoordinate& P);
};

#endif /* ECPRIME_H_ */
