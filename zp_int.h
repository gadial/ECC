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
    zp_int(const mpz_class& _val = 0, const mpz_class& _p = 0):val(_val),p(_p){};
    zp_int(int _val):val(_val),p(0){};
    bool is_equal(const zp_int&) const;
    zp_int& operator=(const zp_int&);
    zp_int& operator+=(const zp_int&);
    zp_int& operator-=(const zp_int&);
    zp_int& operator*=(const zp_int&);
    zp_int& operator/=(const zp_int&);
    zp_int& operator^=(const mpz_class&); //miuse of ^ - means exponent here. WARNING: ^ has low precedence, always use parentheis!
    zp_int& invert();
    zp_int inverse() const;
    ostream& print(ostream&) const;
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

ostream& operator<<(ostream&, const zp_int&);



#endif	/* _ZP_INT_H */

