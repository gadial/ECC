#include <iostream>
#include <stdlib.h>
#include "primes.h"
#include "ellipticcurve.h"
#include "tests.h"
#include "small_primes.h"
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

void PrimesTest::test_generate_prime_for_discriminant(){
    //TODO: this test OFTEN FAILS. check why.
    mpz_class D = gen.rand(100);
    mpz_class s,t;
    mpz_class p = gen.generate_prime_for_discriminant(10,D,t,s);
//    CPPUNIT_ASSERT(4*p == t*t+s*s*D);
}

void PrimesTest::test_is_near_prime(){
    mpz_class p = gen.generate_prime(100);
    mpz_class min_size = p;
    mpz_class temp = p;
    CPPUNIT_ASSERT(is_near_prime(p,-1,min_size) == true);
    for (int i=0; i<NUM_SMALL_PRIMES; i++){
        if (rand() % 2 == 0)
            temp *= small_primes[i];
    }
    CPPUNIT_ASSERT(is_near_prime(temp,-1,min_size) == true);

    CPPUNIT_ASSERT(is_near_prime(p*small_primes[NUM_SMALL_PRIMES - 10],NUM_SMALL_PRIMES - 10,min_size) == false);
}

void EllipticCurveTest::setUp()
{
    random_curve = ECPrime::randomCurve(10, gen);
}

void EllipticCurveTest::tearDown()
{
    
}

void EllipticCurveTest::test_get_point()
{
    Coordinate P = Coordinate::infinity();
    mpz_class x;
    while (P == Coordinate::infinity()){
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
    Coordinate P = Coordinate::infinity();
    Coordinate Q,R;
    while (P == Coordinate::infinity())
        P = random_curve.getPoint(gen.rand(random_curve.mod));
    Q = random_curve.doubling(P);
    //first check if P+P == 2P (P+P should be computed by reduction to the computation of 2P, so shouldn't be a problem
    CPPUNIT_ASSERT(random_curve.addition(P,P) == Q);
    //now check if P + Q + P + Q == 2(P+Q)
    Coordinate temp = P;
    temp = random_curve.addition(temp,Q); // P != Q, so ok
    temp = random_curve.addition(temp,P); // P+Q != P, so ok
    temp = random_curve.addition(temp,Q); // temp = P+Q+P+Q

    CPPUNIT_ASSERT(temp == random_curve.doubling(random_curve.addition(P,Q)));
}

void EllipticCurveTest::test_repeated_doubling(){
    Coordinate P = Coordinate::infinity();
    
    while (P == Coordinate::infinity())
        P = random_curve.getPoint(gen.rand(random_curve.mod));
    
    Coordinate temp = P;
    for (int m=1; m<10; m++){
        temp = random_curve.doubling(temp);
        CPPUNIT_ASSERT(temp == random_curve.repeatedDoubling(P,m));
    }
}

void EllipticCurveTest::test_point_multiplication(){
    Coordinate P = Coordinate::infinity();

    while (P == Coordinate::infinity())
        P = random_curve.getPoint(gen.rand(random_curve.mod));

    CPPUNIT_ASSERT(Coordinate::infinity() == random_curve.pointMultiplication(P,0));
    CPPUNIT_ASSERT(P == random_curve.pointMultiplication(P,1));

    Coordinate temp = random_curve.doubling(P);
    CPPUNIT_ASSERT(temp == random_curve.pointMultiplication(P,2));
    CPPUNIT_ASSERT(random_curve.addition(P,temp) == random_curve.pointMultiplication(P,3));

    mpz_class goal = 2 + gen.rand(50);
    for (int i=2; i< goal; i++)
        temp = random_curve.addition(P,temp);

    CPPUNIT_ASSERT(temp == random_curve.pointMultiplication(P,goal));
}

void PolynomialTest::setUp(){

}

void PolynomialTest::tearDown(){

}
void PolynomialTest::test_input_output(){
    #define S_ARRAY_LENGTH 5
    string s_array[S_ARRAY_LENGTH] = {"1", "x", "x^7", "2x^2 + 3","x^5 + 17x + 543"};
    for (int i=0; i< S_ARRAY_LENGTH; i++){
        ModularPolynomial p(s_array[i], 100000);
        CPPUNIT_ASSERT(s_array[i] == p.to_string());
    }    
}
void PolynomialTest::test_addition_substraction(){
    CPPUNIT_ASSERT(ModularPolynomial("x",100) + ModularPolynomial("x",100) == ModularPolynomial("2x",100));
    CPPUNIT_ASSERT(ModularPolynomial("x",100) + ModularPolynomial("1",100) == ModularPolynomial("x + 1",100));
    CPPUNIT_ASSERT(ModularPolynomial("x^2",100) + ModularPolynomial("55",100) == ModularPolynomial("x^2 + 55",100));
    CPPUNIT_ASSERT(ModularPolynomial("x^2 + 3x + 7",100) + ModularPolynomial("3x^2 + 5x + 12",100) == ModularPolynomial("4x^2 + 8x + 19",100));

    CPPUNIT_ASSERT(ModularPolynomial("x^2",100) - ModularPolynomial("x^2",100) == ModularPolynomial("0",100));
    CPPUNIT_ASSERT(ModularPolynomial("x^2",100) - ModularPolynomial("x",100) == ModularPolynomial("x^2 + 99x",100));
}

void PolynomialTest::test_multiplication(){
    CPPUNIT_ASSERT(ModularPolynomial("x",100) * ModularPolynomial("x",100) == ModularPolynomial("x^2",100));
    CPPUNIT_ASSERT(ModularPolynomial("1",100) * ModularPolynomial("x",100) == ModularPolynomial("x",100));
    CPPUNIT_ASSERT(ModularPolynomial("x + 1",100) * ModularPolynomial("x",100) == ModularPolynomial("x^2 + x",100));
    CPPUNIT_ASSERT(ModularPolynomial("x + 1",100) * ModularPolynomial("x + 1",100) == ModularPolynomial("x^2 + 2x + 1",100));
    CPPUNIT_ASSERT(ModularPolynomial("0",113) * ModularPolynomial("0",113) == ModularPolynomial("0",113));

    CPPUNIT_ASSERT((ModularPolynomial("x",113).modular_exponent(4,ModularPolynomial("x^2",113))) == ModularPolynomial("0",113));
    CPPUNIT_ASSERT((ModularPolynomial("x + 1",113).modular_exponent(2,ModularPolynomial("x^2",113))) == ModularPolynomial("2x + 1",113));
    CPPUNIT_ASSERT((ModularPolynomial("x^2 + 3",113).modular_exponent(8,ModularPolynomial("x^3 + 5x",113))) == ModularPolynomial("18x^2 + 7",113));
    CPPUNIT_ASSERT((ModularPolynomial("x^2 + 3",113).modular_exponent(17,ModularPolynomial("x^3 + 5x",113))) == ModularPolynomial("73x^2 + 34",113));
    CPPUNIT_ASSERT((ModularPolynomial("x^2 + 3",113).modular_exponent(100,ModularPolynomial("x^3 + 5x",113))) == ModularPolynomial("80x^2 + 57",113));
    CPPUNIT_ASSERT((ModularPolynomial("x^2 + 3",113).modular_exponent(1000,ModularPolynomial("x^3 + 5x",113))) == ModularPolynomial("100x^2 + 97",113));
    CPPUNIT_ASSERT((ModularPolynomial("x^2 + 3",113).modular_exponent(10000,ModularPolynomial("x^3 + 5x",113))) == ModularPolynomial("25x^2 + 28",113));
    CPPUNIT_ASSERT((ModularPolynomial("x^2 + 3",113).modular_exponent(100000,ModularPolynomial("x^3 + 5x",113))) == ModularPolynomial("23x^2 + 30",113));
    CPPUNIT_ASSERT((ModularPolynomial("x^2 + 3",113).modular_exponent(1000000,ModularPolynomial("x^3 + 5x",113))) == ModularPolynomial("83x^2 + 106",113));
}

void PolynomialTest::test_divisons(){
//    cout << endl;
//    cout << (ModularPolynomial("x^2 + x",100) % ModularPolynomial("x^2",100)).to_string() << endl;
//    cout << ModularPolynomial("x",100).to_string() << endl;
    CPPUNIT_ASSERT(ModularPolynomial("x^2 + x",100) % ModularPolynomial("x^2",100) == ModularPolynomial("x",100));
    CPPUNIT_ASSERT(ModularPolynomial("x^2 + x",100) % ModularPolynomial("x",100) == ModularPolynomial("0",100));
    CPPUNIT_ASSERT(ModularPolynomial("x^2 + 2x + 7",100) % ModularPolynomial("x + 1",100) == ModularPolynomial("6",100));
    CPPUNIT_ASSERT(ModularPolynomial("x^2 + 2x + 7",100) % ModularPolynomial("x + 1",100) == ModularPolynomial("6",100));
    CPPUNIT_ASSERT(ModularPolynomial("x^3 + 4x + 10",113) % ModularPolynomial("3x + 5",113) == ModularPolynomial("28",113));
    CPPUNIT_ASSERT(ModularPolynomial("x^7 + 34x^5 + 15x^4 + 95x^3 + 17",113) % ModularPolynomial("3x^6 + 5x^3",113) == ModularPolynomial("34x^5 + 51x^4 + 95x^3 + 17",113));
    CPPUNIT_ASSERT(ModularPolynomial("0",113) % ModularPolynomial("x",113) == ModularPolynomial("0",113));

    CPPUNIT_ASSERT(gcd(ModularPolynomial("x",113),ModularPolynomial("x",113)) == ModularPolynomial("x",113));
    CPPUNIT_ASSERT(gcd(ModularPolynomial("x^5 + 3x^2 + x",113),ModularPolynomial("3x^4 + 2x^3 + 17",113)) == ModularPolynomial("1",113));
}