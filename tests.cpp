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
    int number_of_primes = 3;
    int primes[3] = {587, 653, 1033}; // chosen to cover all cases
    for (int i = 0; i<number_of_primes; i++){
        mpz_class n = gen.rand(primes[i]);
        mpz_class root = modular_square_root((n*n) % primes[i],primes[i]);
//        cout << "n, root, p=" << n <<", " << root << ", "<<primes[i]<<endl;
        CPPUNIT_ASSERT(n == root || primes[i]-n == root);
    }
}

void PrimesTest::test_small_primes(){
    for (int i=0; i<10; i++){
        mpz_class p = gen.generate_prime(10);
        for (int j=2; j<p / 2 ; j++) //as naive as it gets
            CPPUNIT_ASSERT(p % j != 0);
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
    random_curve = Ellipticcurve::randomCurve(10, gen);
}

void EllipticCurveTest::tearDown()
{
    
}

void EllipticCurveTest::test_get_point()
{
    Coordinate P = Ellipticcurve::infinity();
    mpz_class x;
    while (P == Ellipticcurve::infinity()){
        x = gen.rand(random_curve.mod);
        P = random_curve.getPoint(x);
    }
    //test if P is really on the curve
//    cout << "x, y = " << P.X << ", " << P.Y << endl;
//    cout << "p, a, b = " << random_curve.mod << ", " << random_curve.ECC_a << ", " <<random_curve.ECC_b << endl;
    CPPUNIT_ASSERT((P.Y*P.Y) % random_curve.mod == (P.X*P.X*P.X + random_curve.ECC_a*P.X + random_curve.ECC_b) % random_curve.mod);
    //test if getPoint really knows to return the "negative" value
    CPPUNIT_ASSERT(P.Y + random_curve.getPoint(x,true).Y == random_curve.mod);
}

void EllipticCurveTest::test_doubling_vs_addition()
{
    Coordinate P = Ellipticcurve::infinity();
    while (P == Ellipticcurve::infinity())
        P = random_curve.getPoint(gen.rand(random_curve.mod));
    //CPPUNIT_ASSERT(random_curve.addition(P,P) == random_curve.doubling(P));
    //TODO: must implement jacobian equality operator in order to test this
}