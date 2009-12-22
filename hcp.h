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
    string to_string();
private:
    map<int,mpz_class> coefficients;
    mpz_class modulos;
    int degree;
};

class HCP{ //Hilbert class polynomials
public:
    HCP();
private:
    map<int, string> H;
};


#endif	/* _HCP_H */

