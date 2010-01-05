/* 
 * File:   main.cpp
 * Author: gadial
 *
 * Created on November 5, 2009, 10:44 AM
 */

#include <iostream>
#include "coordinates.h"
#include "ellipticcurve.h"
#include "curvesnist.h"
#include "primes.h"
#include "curvesnisttest.h"
#include "tests.h"

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
//    cout << (ModularPolynomial("x^3 + 4x + 10",113) % ModularPolynomial("3x + 5",113)) << endl;
//    CPPUNIT_ASSERT() == ModularPolynomial("18x^2 + 7",113));
    do_tests();
    return 0;
}

int do_tests(){

  CPPUNIT_TEST_SUITE_REGISTRATION(CurvesNISTTest);
  CPPUNIT_TEST_SUITE_REGISTRATION(PrimesTest);
  CPPUNIT_TEST_SUITE_REGISTRATION(EllipticCurveTest);
  CPPUNIT_TEST_SUITE_REGISTRATION(PolynomialTest);
  CPPUNIT_TEST_SUITE_REGISTRATION(ZpIntTest);

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
