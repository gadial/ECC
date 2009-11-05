/* 
 * File:   main.cpp
 * Author: gadial
 *
 * Created on November 5, 2009, 10:44 AM
 */

#include <iostream>
#include "coordinates.h"
#include "primes.h"
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    cout << "Hello world!!" << endl;
    cout << "A prime number: " << generate_prime(100) << endl;


    Jacobian jac(1,2,3);
    cout << "X=" << jac.X << ", Y=" << jac.Y << ", Z=" << jac.Z << endl;
    return 0;
}

