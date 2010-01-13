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
    ostream& print(ostream&) const;
    ostream& full_print(ostream&) const;
    operator mpz_class(){return val;}
    mpz_class get_p(){return p;};
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
	ZpCoordinate(mpz_class _x, mpz_class _y, mpz_class _p):
		X(_x,_p), Y(_y,_p),p(_p) {}
	ZpCoordinate(const ZpJacobian& jac);

        //returns the point at infinity, as is represented by this class in the context of elliptic curves
        static ZpCoordinate infinity(){return ZpCoordinate(0,0,0);}

	bool operator==(const ZpCoordinate& eqTo) {
		return X == eqTo.X && Y == eqTo.Y;
	}
        bool isInfinite() {
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
	ZpJacobian(const ZpCoordinate& rhs):
                X(rhs.X), Y(rhs.Y), Z(1,rhs.p) {}

	static ZpJacobian infinity(){return ZpJacobian(1, 1, 0, 0);}
    bool isInfinite() {
    	return X == 1 && Y == 1 && Z == 0;
    }
        mpz_class p;
	zp_int X, Y, Z;
};


#endif	/* _ZP_INT_H */

