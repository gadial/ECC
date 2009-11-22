/* 
 * File:   tests.h
 * Author: gadial
 *
 * Created on November 15, 2009, 10:59 AM
 */

#ifndef _TESTS_H
#define	_TESTS_H

#include <cppunit/extensions/HelperMacros.h>

#include "primes.h"
#include "ellipticcurve.h"
#include "ellipticcurve.h"

class PrimesTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE( PrimesTest );
    CPPUNIT_TEST( test_is_odd );
    CPPUNIT_TEST( test_small_primes );
    CPPUNIT_TEST( test_legendre_symbol );
    CPPUNIT_TEST( test_rand );
    CPPUNIT_TEST( test_square_root );
    CPPUNIT_TEST_SUITE_END();
private:
    RandomNumberGenerator gen;
    mpz_class p,q;
public:
    void setUp();
    void tearDown();

    void test_is_odd();
    void test_legendre_symbol();
    void test_rand();
    void test_square_root();
    void test_small_primes();
};

class EllipticCurveTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE( EllipticCurveTest );
    CPPUNIT_TEST( test_doubling_vs_addition );
    CPPUNIT_TEST( test_get_point );
    CPPUNIT_TEST_SUITE_END();
private:
    RandomNumberGenerator gen;
    Ellipticcurve random_curve;
public:
        void setUp();
        void tearDown();

        void test_get_point();
        void test_doubling_vs_addition();

};

#endif	/* _TESTS_H */

