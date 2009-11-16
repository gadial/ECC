/*
 * curvenisttest.cpp
 *
 *  Created on: Nov 16, 2009
 *      Author: bhess
 */

#include "curvesnisttest.h"

CurvesNISTTest::CurvesNISTTest() {

}

void CurvesNISTTest::setUp() {
	curveP192 = new CurveNISTp192();
	p192S = Coordinate(P192_POINT_S_X, NIST_TESTDATA_BASE,
			P192_POINT_S_Y, NIST_TESTDATA_BASE);
	p192T = Coordinate(P192_POINT_T_X, NIST_TESTDATA_BASE,
			P192_POINT_T_Y, NIST_TESTDATA_BASE);
	p192d.set_str(P192_SCALAR_D, NIST_TESTDATA_BASE);
}

void CurvesNISTTest::tearDown() {
	delete curveP192;
}

void CurvesNISTTest::p192Addition() {
	Coordinate addTestData = Coordinate(P192_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P192_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = Coordinate(curveP192->addition(Jacobian(p192S), p192T),
			curveP192->mod);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p192Subtraction() {
	Coordinate subTestData = Coordinate(P192_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P192_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = Coordinate(curveP192->subtraction(Jacobian(p192S), p192T),
			curveP192->mod);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p192Doubling() {
	Coordinate doublingTestData = Coordinate(P192_2S_X, NIST_TESTDATA_BASE,
			P192_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = Coordinate(curveP192->doubling(Jacobian(p192S)),
			curveP192->mod);
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p192Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P192_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P192_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = Coordinate(curveP192->pointMultiplication(p192S, p192d),
			curveP192->mod);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}
