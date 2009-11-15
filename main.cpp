/* 
 * File:   main.cpp
 * Author: gadial
 *
 * Created on November 5, 2009, 10:44 AM
 */

#include <iostream>
#include "coordinates.h"
#include "primes.h"

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

int do_tests();

using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {

    cout << "Hello world!!" << endl;
    cout << "A prime number: " << generate_prime(100) << endl;


    Jacobian jac(1,2,3);
    cout << "X=" << jac.X << ", Y=" << jac.Y << ", Z=" << jac.Z << endl;

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