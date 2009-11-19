#include <iostream>
#include "primes.h"
#include "ellipticcurve.h"
#include "tests.h"
using namespace std;

void PrimesTest::setUp()
{
    p = gen.generate_prime(100);
    q = 3;
}
void PrimesTest::tearDown()
{
}

void PrimesTest::test_legendre_symbol(){
    int p1 = 587;
    int p1_residues[] = {0, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, 1, 1, -1, 1, -1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, -1, 1, -1, -1, -1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, -1, -1, -1, 1, 1, -1, 1, 1, -1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, 1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, 1, 1, 1, 1, -1, -1, 1, -1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, -1, -1, -1, 1, 1, -1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, -1, 1, -1, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, -1, 1, -1, -1, -1, -1, 1, 1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, -1}; //computed directly via an external script
    for (int i=0; i<p1; i++){
        CPPUNIT_ASSERT(legendre_symbol(i,p1) == p1_residues[i]);
    }

    mpz_class temp_number = gen.rand(p);
    CPPUNIT_ASSERT((temp_number % p == 0) || legendre_symbol(temp_number*temp_number,p) == 1);
    CPPUNIT_ASSERT(legendre_symbol(temp_number*temp_number*temp_number,p) == legendre_symbol(temp_number,p));
    mpz_class mod = p % 3;
    switch(mod.get_ui()){ // -3 is QR modulo p iff p = 1 (mod 3)
        case 0:
        case 2: CPPUNIT_ASSERT(legendre_symbol(p-3,p) == -1); break;
        case 1: CPPUNIT_ASSERT(legendre_symbol(p-3,p) == 1); break;
    }
}

void PrimesTest::test_square_root(){
    int number_of_primes = 2;
    int primes[2] = {587, 653}; // chosen to cover all cases
    for (int i = 0; i<number_of_primes; i++){
        mpz_class n = gen.rand(primes[i]);
        mpz_class root = modular_square_root(n*n,primes[i]);
        cout << "n=" << n <<", root=" << root << " for p="<<primes[i]<<endl;
        CPPUNIT_ASSERT(n == root || primes[i]-n == root);
    }
}

void PrimesTest::test_rand(){
    mpz_class a = gen.rand_binary_digits(100);
    mpz_class b = gen.rand_binary_digits(100);

    CPPUNIT_ASSERT(a != b); //if this happens, our PRNG sucks.
}

void PrimesTest::test_is_odd()
{
    CPPUNIT_ASSERT((p % 2) == 1);
    CPPUNIT_ASSERT((q % 2) == 1);
}

void EllipticCurveTest::setUp()
{

}

void EllipticCurveTest::tearDown()
{
    
}
void EllipticCurveTest::test_doubling_vs_addition()
{

}