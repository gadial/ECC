/*
 * padictest.h
 *
 *  Created on: Jan 23, 2010
 *      Author: bhess
 */

#ifndef PADICTEST_H_
#define PADICTEST_H_

#include <cppunit/extensions/HelperMacros.h>
#include "../arith/Poly.h"
#include "../adicops.h"

class Padictest : public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE(Padictest);
	CPPUNIT_TEST(teichmueller_mod);
	CPPUNIT_TEST(fast_division_with_remainder);
	CPPUNIT_TEST(inverse);
	CPPUNIT_TEST(invsqrt);
	CPPUNIT_TEST_SUITE_END();
public:
	Padictest() {};
	void setUp();
	void tearDown();

private:

	/*
	 * According to Example 12.4
	 */
	void teichmueller_mod();

	/*
	 * According to Example 12.7
	 */
	void fast_division_with_remainder();

	/*
	 * According to Example 12.11
	 */
	void inverse();

	/*
	 * According to Example 12.13
	 */
	void invsqrt();

	Adicops* ad;
	Poly* M;
};

#endif /* PADICTEST_H_ */
