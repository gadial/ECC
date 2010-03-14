/*
 * bincurvetest.cpp
 *
 *  Created on: Mar 10, 2010
 *      Author: bhess
 */

#include "bincurvetest.h"

void BinCurveTest::setUp() {
	ecb = new CurveNISTb163();
}

void BinCurveTest::tearDown() {
	delete ecb;
}

void BinCurveTest::compressed_format_check() {
	// Compressing and uncompressing a point...
	//cout << ecb->point << endl;
	string comp = ecb->toCompressedForm(ecb->point);
	//cout << comp << endl;
	Coordinate co = ecb->getPointCompressedForm(comp);
	//cout << co << endl;
	//ecb->check_coordinate(co);
	CPPUNIT_ASSERT(ecb->point == co);
}

void BinCurveTest::check_order() {
	Coordinate m = ecb->pointMultiplication(ecb->point, ecb->getOrder());
	CPPUNIT_ASSERT(m.isInfinite());
}
