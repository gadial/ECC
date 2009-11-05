/*
 * ellipticcurve.cpp
 *
 *  Created on: Nov 5, 2009
 *      Author: bhess
 */

#include "ellipticcurve.h"

Ellipticcurve::Ellipticcurve() {}

Ellipticcurve::~Ellipticcurve() {}

Jacobian Ellipticcurve::addition(Jacobian P, Coordinate Q) {}

Jacobian Ellipticcurve::doubling(Jacobian P) {}

Coordinate Ellipticcurve::pointMultiplication(Coordinate P, mpz_class k) {}

Jacobian Ellipticcurve::repeatedDoubling(Jacobian P, int m) {}

