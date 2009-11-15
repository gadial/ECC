/* 
 * File:   main.cpp
 * Author: gadial
 *
 * Created on November 5, 2009, 10:44 AM
 */

#include <iostream>
#include "coordinates.h"
#include "ellipticcurve.h"
#include "primes.h"

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

int do_tests();

using namespace std;

void printCd(const Coordinate& c) {
	cout << "(" << c.X << "," << c.Y << ")" << endl;
}
void printJac(const Jacobian& j) {
	cout << "(" << j.X << "," << j.Y << "," << j.Z << ")" << endl;
}

/*
 * 
 */
int main(int argc, char** argv) {

    cout << "Hello world!!" << endl;
    cout << "A prime number: " << generate_prime(100) << endl;

    /*
    Coordinate c1(27,27);
    Jacobian j1(c1);
    Coordinate c2(j1, 29);
    printCd(c1); printJac(j1); printCd(c1);
    */

    Ellipticcurve* ellC = new Ellipticcurve(29, 4, 20);
    Jacobian add = ellC->addition(Jacobian(1, 1, 0), Coordinate(16, 27));
    Jacobian doub = ellC->doubling(Jacobian(Coordinate(1,5)));
    Jacobian mult = ellC->pointMultiplication(Coordinate(1,5), 33);


    printCd(Coordinate(add, 29));
    printCd(Coordinate(doub, 29));
    printCd(Coordinate(mult, 29));

    cout << "X=" << add.X << ", Y=" << add.Y << ", Z=" << add.Z << endl;

    do_tests();
    return 0;
}

int do_tests(){
    // Get the top level suite from the registry
  CppUnit::Test *suite = CppUnit::TestFactoryRegistry::getRegistry().makeTest();

  // Adds the test to the list of test to run
  CppUnit::TextUi::TestRunner runner;
  runner.addTest( suite );

  // Change the default outputter to a compiler error format outputter
  runner.setOutputter( new CppUnit::CompilerOutputter( &runner.result(),
                                                       std::cerr ) );
  // Run the tests.
  bool wasSucessful = runner.run();

  // Return error code 1 if the one of test failed.
  return wasSucessful ? 0 : 1;
}