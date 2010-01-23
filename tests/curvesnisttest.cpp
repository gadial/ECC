/*
 * curvenisttest.cpp
 *
 *  Created on: Nov 16, 2009
 *      Author: bhess
 */

#include "curvesnisttest.h"
#include "../zp_int.h"

CurvesNISTTest::CurvesNISTTest() {

}

void CurvesNISTTest::setUp() {
	curveP192 = new CurveNISTp192();
	p192S = Coordinate(P192_POINT_S_X, NIST_TESTDATA_BASE,
			P192_POINT_S_Y, NIST_TESTDATA_BASE);
	p192T = Coordinate(P192_POINT_T_X, NIST_TESTDATA_BASE,
			P192_POINT_T_Y, NIST_TESTDATA_BASE);
	p192d.set_str(P192_SCALAR_D, NIST_TESTDATA_BASE);

	curveP224 = new CurveNISTp224();
	p224S = Coordinate(P224_POINT_S_X, NIST_TESTDATA_BASE,
			P224_POINT_S_Y, NIST_TESTDATA_BASE);
	p224T = Coordinate(P224_POINT_T_X, NIST_TESTDATA_BASE,
			P224_POINT_T_Y, NIST_TESTDATA_BASE);
	p224d.set_str(P224_SCALAR_D, NIST_TESTDATA_BASE);

	curveP256 = new CurveNISTp256();
	p256S = Coordinate(P256_POINT_S_X, NIST_TESTDATA_BASE,
			P256_POINT_S_Y, NIST_TESTDATA_BASE);
	p256T = Coordinate(P256_POINT_T_X, NIST_TESTDATA_BASE,
			P256_POINT_T_Y, NIST_TESTDATA_BASE);
	p256d.set_str(P256_SCALAR_D, NIST_TESTDATA_BASE);

	curveP384 = new CurveNISTp384();
	p384S = Coordinate(P384_POINT_S_X, NIST_TESTDATA_BASE,
			P384_POINT_S_Y, NIST_TESTDATA_BASE);
	p384T = Coordinate(P384_POINT_T_X, NIST_TESTDATA_BASE,
			P384_POINT_T_Y, NIST_TESTDATA_BASE);
	p384d.set_str(P384_SCALAR_D, NIST_TESTDATA_BASE);

	curveP521 = new CurveNISTp521();
	p521S = Coordinate(P521_POINT_S_X, NIST_TESTDATA_BASE,
			P521_POINT_S_Y, NIST_TESTDATA_BASE);
	p521T = Coordinate(P521_POINT_T_X, NIST_TESTDATA_BASE,
			P521_POINT_T_Y, NIST_TESTDATA_BASE);
	p521d.set_str(P521_SCALAR_D, NIST_TESTDATA_BASE);

	curveB163 = new CurveNISTb163();
}

void CurvesNISTTest::tearDown() {
	delete curveP192;
	delete curveP224;
	delete curveP256;
	delete curveP384;
	delete curveP521;
	delete curveB163;
}

void CurvesNISTTest::p192Addition() {
	Coordinate addTestData = Coordinate(P192_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P192_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = curveP192->addition(p192S, p192T);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p192Subtraction() {
	Coordinate subTestData = Coordinate(P192_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P192_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = curveP192->subtraction(p192S, p192T);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p192Doubling() {
	Coordinate doublingTestData = Coordinate(P192_2S_X, NIST_TESTDATA_BASE,
			P192_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = curveP192->doubling(p192S);
        zp_int x1(p192S.X,curveP192->mod);
        zp_int y1(p192S.Y,curveP192->mod);
        zp_int x3 = ((((x1^2)*3+curveP192->ECC_a)^2)-x1*(y1^2)*8)/((y1^2)*4);
        zp_int y3 = (((((x1^2)*3+curveP192->ECC_a))/(y1*2))*(x1-x3))-y1;       
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p192Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P192_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P192_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = curveP192->pointMultiplication(p192S, p192d);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::p192Order() {
	Coordinate multiplicationCalc = curveP192->pointMultiplication(curveP192->point, curveP192->getOrder());
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}

void CurvesNISTTest::p224Addition() {
	Coordinate addTestData = Coordinate(P224_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P224_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = curveP224->addition(p224S, p224T);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p224Subtraction() {
	Coordinate subTestData = Coordinate(P224_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P224_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = curveP224->subtraction(p224S, p224T);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p224Doubling() {
	Coordinate doublingTestData = Coordinate(P224_2S_X, NIST_TESTDATA_BASE,
			P224_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = curveP224->doubling(p224S);
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p224Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P224_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P224_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = curveP224->pointMultiplication(p224S, p224d);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::p224Order() {
	Coordinate multiplicationCalc = curveP224->pointMultiplication(curveP224->point, curveP224->getOrder());
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}

void CurvesNISTTest::p256Addition() {
	Coordinate addTestData = Coordinate(P256_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P256_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = curveP256->addition(p256S, p256T);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p256Subtraction() {
	Coordinate subTestData = Coordinate(P256_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P256_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = curveP256->subtraction(p256S, p256T);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p256Doubling() {
	Coordinate doublingTestData = Coordinate(P256_2S_X, NIST_TESTDATA_BASE,
			P256_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = curveP256->doubling(p256S);

	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p256Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P256_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P256_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = curveP256->pointMultiplication(p256S, p256d);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::p256Order() {
	Coordinate multiplicationCalc = curveP256->pointMultiplication(curveP256->point, curveP256->getOrder());
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}

void CurvesNISTTest::p384Addition() {
	Coordinate addTestData = Coordinate(P384_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P384_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = curveP384->addition(p384S, p384T);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p384Subtraction() {
	Coordinate subTestData = Coordinate(P384_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P384_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = curveP384->subtraction(p384S, p384T);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p384Doubling() {
	Coordinate doublingTestData = Coordinate(P384_2S_X, NIST_TESTDATA_BASE,
			P384_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = curveP384->doubling(p384S);
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p384Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P384_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P384_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = curveP384->pointMultiplication(p384S, p384d);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::p384Order() {
	Coordinate multiplicationCalc = curveP384->pointMultiplication(curveP384->point, curveP384->getOrder());
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}

void CurvesNISTTest::p521Addition() {
	Coordinate addTestData = Coordinate(P521_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P521_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = curveP521->addition(p521S, p521T);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p521Subtraction() {
	Coordinate subTestData = Coordinate(P521_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P521_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = curveP521->subtraction(p521S, p521T);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p521Doubling() {
	Coordinate doublingTestData = Coordinate(P521_2S_X, NIST_TESTDATA_BASE,
			P521_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = curveP521->doubling(p521S);
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p521Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P521_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P521_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = curveP521->pointMultiplication(p521S, p521d);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::p521Order() {
	Coordinate multiplicationCalc = curveP521->pointMultiplication(curveP521->point, curveP521->getOrder());
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}

void CurvesNISTTest::b163Addition() {
	// TODO: find testdata
	Coordinate addTestData = Coordinate(P521_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P521_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = curveP521->addition(p521S, p521T);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::b163Subtraction() {
	// TODO: find testdata
	Coordinate subTestData = Coordinate(P521_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P521_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = curveP521->subtraction(p521S, p521T);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::b163Doubling() {
	// TODO: find testdata
	Coordinate doublingTestData = Coordinate(P521_2S_X, NIST_TESTDATA_BASE,
			P521_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = curveP521->doubling(p521S);
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::b163Multiplication() {
	// TODO: find testdata
	Coordinate multiplicationTestData = Coordinate(P521_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P521_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = curveP521->pointMultiplication(p521S, p521d);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::b163Order() {
	//std::cout << curveB163->mod << std::endl;
	Coordinate multiplicationCalc = curveB163->pointMultiplication(curveB163->point, curveB163->getOrder());
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}
