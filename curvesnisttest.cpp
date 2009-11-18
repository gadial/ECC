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
}

void CurvesNISTTest::tearDown() {
	delete curveP192;
	delete curveP224;
	delete curveP256;
	delete curveP384;
	delete curveP521;
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

void CurvesNISTTest::p192Order() {
	Coordinate multiplicationCalc = Coordinate(curveP192->pointMultiplication(curveP192->point, curveP192->getOrder()),
			curveP192->mod);
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}

void CurvesNISTTest::p224Addition() {
	Coordinate addTestData = Coordinate(P224_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P224_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = Coordinate(curveP224->addition(Jacobian(p224S), p224T),
			curveP224->mod);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p224Subtraction() {
	Coordinate subTestData = Coordinate(P224_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P224_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = Coordinate(curveP224->subtraction(Jacobian(p224S), p224T),
			curveP224->mod);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p224Doubling() {
	Coordinate doublingTestData = Coordinate(P224_2S_X, NIST_TESTDATA_BASE,
			P224_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = Coordinate(curveP224->doubling(Jacobian(p224S)),
			curveP224->mod);
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p224Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P224_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P224_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = Coordinate(curveP224->pointMultiplication(p224S, p224d),
			curveP224->mod);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::p224Order() {
	Coordinate multiplicationCalc = Coordinate(curveP224->pointMultiplication(curveP224->point, curveP224->getOrder()),
			curveP224->mod);
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}

void CurvesNISTTest::p256Addition() {
	Coordinate addTestData = Coordinate(P256_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P256_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = Coordinate(curveP256->addition(Jacobian(p256S), p256T),
			curveP256->mod);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p256Subtraction() {
	Coordinate subTestData = Coordinate(P256_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P256_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = Coordinate(curveP256->subtraction(Jacobian(p256S), p256T),
			curveP256->mod);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p256Doubling() {
	Coordinate doublingTestData = Coordinate(P256_2S_X, NIST_TESTDATA_BASE,
			P256_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = Coordinate(curveP256->doubling(Jacobian(p256S)),
			curveP256->mod);
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p256Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P256_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P256_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = Coordinate(curveP256->pointMultiplication(p256S, p256d),
			curveP256->mod);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::p256Order() {
	Coordinate multiplicationCalc = Coordinate(curveP256->pointMultiplication(curveP256->point, curveP256->getOrder()),
			curveP256->mod);
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}

void CurvesNISTTest::p384Addition() {
	Coordinate addTestData = Coordinate(P384_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P384_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = Coordinate(curveP384->addition(Jacobian(p384S), p384T),
			curveP384->mod);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p384Subtraction() {
	Coordinate subTestData = Coordinate(P384_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P384_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = Coordinate(curveP384->subtraction(Jacobian(p384S), p384T),
			curveP384->mod);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p384Doubling() {
	Coordinate doublingTestData = Coordinate(P384_2S_X, NIST_TESTDATA_BASE,
			P384_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = Coordinate(curveP384->doubling(Jacobian(p384S)),
			curveP384->mod);
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p384Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P384_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P384_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = Coordinate(curveP384->pointMultiplication(p384S, p384d),
			curveP384->mod);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::p384Order() {
	Coordinate multiplicationCalc = Coordinate(curveP384->pointMultiplication(curveP384->point, curveP384->getOrder()),
			curveP384->mod);
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}

void CurvesNISTTest::p521Addition() {
	Coordinate addTestData = Coordinate(P521_S_PLUS_T_X, NIST_TESTDATA_BASE,
			P521_S_PLUS_T_Y, NIST_TESTDATA_BASE);
	Coordinate addCalc = Coordinate(curveP521->addition(Jacobian(p521S), p521T),
			curveP521->mod);

	CPPUNIT_ASSERT(addTestData == addCalc);
}

void CurvesNISTTest::p521Subtraction() {
	Coordinate subTestData = Coordinate(P521_S_MINUS_T_X, NIST_TESTDATA_BASE,
			P521_S_MINUS_T_Y, NIST_TESTDATA_BASE);

	Coordinate subCalc = Coordinate(curveP521->subtraction(Jacobian(p521S), p521T),
			curveP521->mod);

	CPPUNIT_ASSERT(subTestData == subCalc);
}

void CurvesNISTTest::p521Doubling() {
	Coordinate doublingTestData = Coordinate(P521_2S_X, NIST_TESTDATA_BASE,
			P521_2S_Y, NIST_TESTDATA_BASE);
	Coordinate doublingCalc = Coordinate(curveP521->doubling(Jacobian(p521S)),
			curveP521->mod);
	CPPUNIT_ASSERT(doublingTestData == doublingCalc);
}

void CurvesNISTTest::p521Multiplication() {
	Coordinate multiplicationTestData = Coordinate(P521_D_TIMES_S_X, NIST_TESTDATA_BASE,
			P521_D_TIMES_S_Y, NIST_TESTDATA_BASE);
	Coordinate multiplicationCalc = Coordinate(curveP521->pointMultiplication(p521S, p521d),
			curveP521->mod);
	CPPUNIT_ASSERT(multiplicationTestData == multiplicationCalc);
}

void CurvesNISTTest::p521Order() {
	Coordinate multiplicationCalc = Coordinate(curveP521->pointMultiplication(curveP521->point, curveP521->getOrder()),
			curveP521->mod);
	CPPUNIT_ASSERT(multiplicationCalc.isInfinite());
}
