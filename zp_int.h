/* 
 * File:   zp_int.h
 * Author: gadial
 *
 * Created on January 3, 2010, 3:41 PM
 */

//zp_int is an element of Z_p, for some prime p
//this class wraps the basic arithmetic functions of such numbers, in order to simplify algorithms using it.
//We rely on the gmp library in order to allow arbitarly-large numbers
#ifndef _ZP_INT_H
#define	_ZP_INT_H
#include <iostream>
#include <gmpxx.h>
#include "coordinates.h"
using std::ostream;

class zp_int{
public:
    zp_int(const mpz_class& _val = 0, const mpz_class& _p = 0):val(_val),p(_p){normalize();};
    zp_int(int _val):val(_val),p(0){};
    bool is_equal(const zp_int&) const;
    bool is_smaller(const zp_int&) const;
    zp_int& operator=(const zp_int&);
    zp_int& operator+=(const zp_int&);
    zp_int& operator-=(const zp_int&);
    zp_int operator-() const;
    zp_int& operator*=(const zp_int&);
    zp_int& operator/=(const zp_int&);
    zp_int& operator^=(const mpz_class&); //miuse of ^ - means exponent here. WARNING: ^ has low precedence, always use parentheis!
    zp_int& invert();
    zp_int inverse() const;
    int operator%(const int n) const {mpz_class temp = val % n; return temp.get_ui();}
    ostream& print(ostream&) const;
    ostream& full_print(ostream&) const;
    operator mpz_class(){return val;}
    mpz_class get_p() const {return p;};
    string to_s(int base = 10) const;
private:
    mpz_class val;
    mpz_class p; // p=0 means we treat it as a normal integer

    zp_int& normalize(); //sets val to be the equivalent value mod p
};

zp_int operator+(const zp_int&, const zp_int&);
zp_int operator-(const zp_int&, const zp_int&);
zp_int operator*(const zp_int&, const zp_int&);
zp_int operator/(const zp_int&, const zp_int&);
zp_int operator^(const zp_int&, const mpz_class&);
bool operator==(const zp_int&, const zp_int&);
bool operator!=(const zp_int&, const zp_int&);
bool operator<(const zp_int&, const zp_int&);

ostream& operator<<(ostream&, const zp_int&);


//zp_int based coordinates
class ZpCoordinate;
class ZpJacobian;

class ZpCoordinate {
public:

	ZpCoordinate() {}
//	ZpCoordinate(mpz_class _x, mpz_class _y, mpz_class _p):
//		X(_x,_p), Y(_y,_p),p(_p) {}
        ZpCoordinate(zp_int _x, zp_int _y, mpz_class _p):
		X(_x), Y(_y),p(_p) {}
	ZpCoordinate(const ZpJacobian& jac);
        ZpCoordinate(const Coordinate& cor, mpz_class _p):
            X(cor.X,_p),Y(cor.Y,_p),p(_p){}
        operator Coordinate(){return Coordinate(X,Y);}
        string toCompressedForm();
        //returns the point at infinity, as is represented by this class in the context of elliptic curves
        static ZpCoordinate infinity(){return ZpCoordinate(0,0,0);}

	bool operator==(const ZpCoordinate& eqTo) {
		return X == eqTo.X && Y == eqTo.Y;
	}
        bool isInfinite() const {
        	return X == 0 && Y == 0;
        }
        mpz_class p;
	zp_int X, Y;
};

ostream& operator<<(ostream& out, const ZpCoordinate& rhs);

class ZpJacobian {
public:

	ZpJacobian() {}
	ZpJacobian(mpz_class _x, mpz_class _y, mpz_class _z, mpz_class _p):
		X(_x, _p), Y(_y, _p), Z(_z, _p), p(_p) {}
        ZpJacobian(zp_int _x, zp_int _y, zp_int _z):
        	X(_x), Y(_y), Z(_z), p(_x.get_p()) {}
	ZpJacobian(const ZpCoordinate& rhs):
                X(rhs.X), Y(rhs.Y), Z(1,rhs.p),p(rhs.p) {}

	static ZpJacobian infinity(mpz_class p = 0){return ZpJacobian(1, 1, 0, p);}
        operator Jacobian(){return Jacobian(X,Y,Z);}
    bool isInfinite() const {
    	return X == 1 && Y == 1 && Z == 0;
    }
        mpz_class p;
	zp_int X, Y, Z;
};

ostream& operator<<(ostream& out, const ZpJacobian& rhs);

#endif	/* _ZP_INT_H */

