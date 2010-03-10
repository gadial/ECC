/*
 * bincurvetest.cpp
 *
 *  Created on: Mar 10, 2010
 *      Author: bhess
 */

#include "bincurvetest.h"
#include "../curvesnist.h"

void BinCurveTest::setUp() {

}

void BinCurveTest::tearDown() {

}

void BinCurveTest::compressed_format_check() {
	ECBinary* ecb = new CurveNISTb163();

	cout << ecb->point << endl;
	string comp = ecb->toCompressedForm(ecb->point);
	cout << comp << endl;
	Coordinate co = ecb->getPointCompressedForm(comp);
	cout << co << endl;
	//ecb->check_coordinate(co);

	delete ecb;
}
