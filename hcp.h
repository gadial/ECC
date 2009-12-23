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

class ModularPolynomial{
public:
    ModularPolynomial(string, mpz_class);
    ModularPolynomial(const ModularPolynomial&);
    string to_string();
    ModularPolynomial& operator=(const ModularPolynomial&);
    ModularPolynomial& operator+=(const ModularPolynomial&);
    ModularPolynomial& operator-=(const ModularPolynomial&);
    ModularPolynomial& operator*=(const ModularPolynomial&);
    bool operator==(const ModularPolynomial&) const;
private:
    mutable map<int,mpz_class> coefficients; //mutable since otherwise we get const-correctness problem (seems like mpz_class is not very const-correct)
    mpz_class modulus;
    int degree;
};

ModularPolynomial operator+(const ModularPolynomial& rhs, const ModularPolynomial& lhs);
ModularPolynomial operator-(const ModularPolynomial& rhs, const ModularPolynomial& lhs);
ModularPolynomial operator*(const ModularPolynomial& rhs, const ModularPolynomial& lhs);

class HCP{ //Hilbert class polynomials
public:
    HCP();
private:
    map<int, string> H;
};


#endif	/* _HCP_H */

