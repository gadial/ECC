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
#include "tests/curvesnisttest.h"
#include "tests/padictest.h"
#include "tests.h"
#include "adicops.h"
#include "elgamal.h"
#include "cmd.h"
#include "challange_crack.h"

//#include <NTL/GF2EX.h>

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
    try_all_curves_on_cipher("bla", "bla", "bla");
//    	Cmd* cmd = new Cmd(argc, argv);
//	if (cmd->do_tests) {
//		do_tests();
//	} else {
//		cmd->execute();
//	}
//	delete cmd;
    return 0;
}

int do_tests(){

  CPPUNIT_TEST_SUITE_REGISTRATION(CurvesNISTTest);
  CPPUNIT_TEST_SUITE_REGISTRATION(PrimesTest);
  CPPUNIT_TEST_SUITE_REGISTRATION(EllipticCurveTest);
  CPPUNIT_TEST_SUITE_REGISTRATION(PolynomialTest);
  CPPUNIT_TEST_SUITE_REGISTRATION(ZpIntTest);
  CPPUNIT_TEST_SUITE_REGISTRATION(Padictest);

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
