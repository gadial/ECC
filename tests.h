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
#include "ecprime.h"
#include "hcp.h"
#include "zp_int.h"

class PrimesTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE( PrimesTest );
    CPPUNIT_TEST( test_is_odd );
    CPPUNIT_TEST( test_small_primes );
    CPPUNIT_TEST( test_legendre_symbol );
    CPPUNIT_TEST( test_rand );
    CPPUNIT_TEST( test_square_root );
    CPPUNIT_TEST( test_generate_prime_for_discriminant );
    CPPUNIT_TEST( test_is_near_prime );
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
    void test_generate_prime_for_discriminant();
    void test_is_near_prime();
};

class EllipticCurveTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE( EllipticCurveTest );
    CPPUNIT_TEST( test_doubling_vs_addition );
    CPPUNIT_TEST( test_get_point );
    CPPUNIT_TEST( test_repeated_doubling );
    CPPUNIT_TEST( test_point_multiplication );
    CPPUNIT_TEST( test_check_order );
    CPPUNIT_TEST( test_coordinate_compressed_form );
    CPPUNIT_TEST_SUITE_END();
private:
    RandomNumberGenerator gen;
    ECPrime random_curve;
public:
        void setUp();
        void tearDown();

        void test_get_point();
        void test_doubling_vs_addition();
        void test_repeated_doubling();
        void test_point_multiplication();
        void test_check_order();
        void test_coordinate_compressed_form();

};

class PolynomialTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE( PolynomialTest );
//    CPPUNIT_TEST( test_input_output );
//    CPPUNIT_TEST( test_addition_substraction );
//    CPPUNIT_TEST( test_multiplication );
//    CPPUNIT_TEST( test_divisons );
//    CPPUNIT_TEST( test_evaluations );
    CPPUNIT_TEST( test_root_finding );
    CPPUNIT_TEST_SUITE_END();
private:
    #define ROOTS_ARRAY_LENGTH 8
    RandomNumberGenerator gen;
    mpz_class p;
    NumberArray random_roots;
public:
    void setUp();
    void tearDown();

    void test_input_output();
    void test_addition_substraction();
    void test_multiplication();
    void test_divisons();
    void test_evaluations();
    void test_root_finding();
};

class ZpIntTest : public CppUnit::TestFixture{
    CPPUNIT_TEST_SUITE( ZpIntTest );
    CPPUNIT_TEST( test_arithmetic );
    CPPUNIT_TEST_SUITE_END();
private:
    #define NUMBER_ARRAY_LENGTH 100
    RandomNumberGenerator gen;
    zp_int numbers[NUMBER_ARRAY_LENGTH];
public:
    void setUp();
    void tearDown();

    void test_arithmetic();
};
#endif	/* _TESTS_H */

