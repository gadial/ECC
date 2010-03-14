/*
 * bincurvetest.h
 *
 *  Created on: Mar 10, 2010
 *      Author: bhess
 */

#ifndef BINCURVETEST_H_
#define BINCURVETEST_H_

#include <cppunit/extensions/HelperMacros.h>
#include "../curvesnist.h"

class BinCurveTest  : public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE(BinCurveTest);
	CPPUNIT_TEST(compressed_format_check);
	CPPUNIT_TEST(check_order);
	CPPUNIT_TEST_SUITE_END();
public:
	BinCurveTest() {};
	void setUp();
	void tearDown();

private:
	ECBinary* ecb;

	/*
	 * According to Example 12.4
	 */
	void compressed_format_check();

	void check_order();

};

#endif /* BINCURVETEST_H_ */
