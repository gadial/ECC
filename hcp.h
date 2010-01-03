/* 
 * File:   hcp.h
 * Author: gadial
 *
 * Created on December 3, 2009, 11:50 AM
 */

#ifndef _HCP_H
#define	_HCP_H

#include <map>
#include <string>
#include <gmpxx.h>
#include <vector>
#include <iostream>

using std::map;
using std::string;
using std::vector;
using std::ostream;

typedef vector<mpz_class> NumberArray;

class ModularPolynomial{
public:
    ModularPolynomial():modulus(0),degree(-1){};
    ModularPolynomial(string, mpz_class);
    ModularPolynomial(const ModularPolynomial&);
    string to_string() const;
    ModularPolynomial& operator=(const ModularPolynomial&);
    ModularPolynomial& operator+=(const ModularPolynomial&);
    ModularPolynomial& operator-=(const ModularPolynomial&);
    ModularPolynomial& operator*=(const ModularPolynomial&);
    ModularPolynomial& operator/=(const ModularPolynomial&);
    ModularPolynomial modular_exponent(mpz_class exp, const ModularPolynomial& mod) const;
    mpz_class operator()(mpz_class a) const;

    ModularPolynomial operator%(const ModularPolynomial& lhs);
    bool operator==(const ModularPolynomial&) const;
    ModularPolynomial& normalize();
    bool is_zero() const;
    NumberArray find_roots();
    int get_degree(){return degree;}
private:
    mutable map<int,mpz_class> coefficients; //mutable since otherwise we get const-correctness problem (seems like mpz_class is not very const-correct)
    mpz_class modulus;
    int degree;

    void divide(ModularPolynomial lhs, ModularPolynomial& Q,ModularPolynomial& R);
};

ModularPolynomial operator+(const ModularPolynomial& rhs, const ModularPolynomial& lhs);
ModularPolynomial operator-(const ModularPolynomial& rhs, const ModularPolynomial& lhs);
ModularPolynomial operator*(const ModularPolynomial& rhs, const ModularPolynomial& lhs);
ModularPolynomial operator/(const ModularPolynomial& rhs, const ModularPolynomial& lhs);
ModularPolynomial gcd(const ModularPolynomial& rhs, const ModularPolynomial& lhs);


ostream& operator<<(ostream& o, const ModularPolynomial& lhs);
ostream& operator<<(ostream& o, const NumberArray lhs);

class HCP{ //Hilbert class polynomials
public:
    HCP();
private:
    map<int, string> H;
};


#endif	/* _HCP_H */

