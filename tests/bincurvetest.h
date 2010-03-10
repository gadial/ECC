/*
 * bincurvetest.h
 *
 *  Created on: Mar 10, 2010
 *      Author: bhess
 */

#ifndef BINCURVETEST_H_
#define BINCURVETEST_H_

#include <cppunit/extensions/HelperMacros.h>

class BinCurveTest  : public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE(BinCurveTest);
	CPPUNIT_TEST(compressed_format_check);
	CPPUNIT_TEST_SUITE_END();
public:
	BinCurveTest() {};
	void setUp();
	void tearDown();

private:

	/*
	 * According to Example 12.4
	 */
	void compressed_format_check();

};

#endif /* BINCURVETEST_H_ */
