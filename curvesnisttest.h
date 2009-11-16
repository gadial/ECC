/*
 * curvesnisttest.h
 *
 *  Created on: Nov 16, 2009
 *      Author: bhess
 *
 *  Testsuite for the NIST FIPS 186-3 ECs
 *
 *  Test data from NSA Suite B
 *  http://www.nsa.gov/ia/_files/nist-routines.pdf
 *
 */

#ifndef CURVESNISTTEST_H_
#define CURVESNISTTEST_H_

#include <cppunit/extensions/HelperMacros.h>
#include "coordinates.h"
#include "curvesnist.h"

#define NIST_TESTDATA_BASE 16

#define P192_POINT_S_X "d458e7d1 27ae671b 0c330266 d2467693 53a01207 3e97acf8"
#define P192_POINT_S_Y "32593050 0d851f33 6bddc050 cf7fb11b 5673a164 5086df3b"
#define P192_POINT_T_X "f22c4395 213e9ebe 67ddecdd 87fdbd01 be16fb05 9b9753a4"
#define P192_POINT_T_Y "26442409 6af2b359 7796db48 f8dfb41f a9cecc97 691a9c79"
#define P192_S_PLUS_T_X "48e1e409 6b9b8e5c a9d0f1f0 77b8abf5 8e843894 de4d0290"
#define P192_S_PLUS_T_Y "408fa77c 797cd7db fb16aa48 a3648d3d 63c94117 d7b6aa4b"
#define P192_S_MINUS_T_X "fc9683cc 5abfb4fe 0cc8cc3b c9f61eab c4688f11 e9f64a2e"
#define P192_S_MINUS_T_Y "093e31d0 0fb78269 732b1bd2 a73c23cd d31745d0 523d816b"
#define P192_2S_X "30c5bc6b 8c7da253 54b373dc 14dd8a0e ba42d25a 3f6e6962"
#define P192_2S_Y "0dde14bc 4249a721 c407aedb f011e2dd bbcb2968 c9d889cf"
#define P192_SCALAR_D "a78a236d 60baec0c 5dd41b33 a542463a 8255391a f64c74ee"
#define P192_D_TIMES_S_X "1faee420 5a4f669d 2d0a8f25 e3bcec9a 62a69529 65bf6d31"
#define P192_D_TIMES_S_Y "5ff2cdfa 508a2581 89236708 7c696f17 9e7a4d7e 8260fb06"


class CurvesNISTTest : public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE(CurvesNISTTest);
	CPPUNIT_TEST(p192Addition);
	CPPUNIT_TEST(p192Subtraction);
	CPPUNIT_TEST(p192Doubling);
	CPPUNIT_TEST(p192Multiplication);
	CPPUNIT_TEST_SUITE_END();

public:
	CurvesNISTTest();

	void setUp();
	void tearDown();

	void p192Addition();
	void p192Subtraction();
	void p192Doubling();
	void p192Multiplication();

private:
	Ellipticcurve* curveP192;
	Coordinate p192S, p192T;
	mpz_class p192d;

};

#endif /* CURVESNISTTEST_H_ */
