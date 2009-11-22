/*
 * curvesnisttest.h
 *
 *  Created on: Nov 16, 2009
 *      Author: bhess
 *
 *  Testsuite for the NIST FIPS 186-3 ECs
 *
 *  Test data from NSA Suite B
 *  http://www.nsa.gov/ia/_files/nist-routines.pdf
 *
 */

#ifndef CURVESNISTTEST_H_
#define CURVESNISTTEST_H_

#include <cppunit/extensions/HelperMacros.h>
#include "curvesnist.h"

#define NIST_TESTDATA_BASE 16

#define P192_POINT_S_X "d458e7d1 27ae671b 0c330266 d2467693 53a01207 3e97acf8"
#define P192_POINT_S_Y "32593050 0d851f33 6bddc050 cf7fb11b 5673a164 5086df3b"
#define P192_POINT_T_X "f22c4395 213e9ebe 67ddecdd 87fdbd01 be16fb05 9b9753a4"
#define P192_POINT_T_Y "26442409 6af2b359 7796db48 f8dfb41f a9cecc97 691a9c79"
#define P192_S_PLUS_T_X "48e1e409 6b9b8e5c a9d0f1f0 77b8abf5 8e843894 de4d0290"
#define P192_S_PLUS_T_Y "408fa77c 797cd7db fb16aa48 a3648d3d 63c94117 d7b6aa4b"
#define P192_S_MINUS_T_X "fc9683cc 5abfb4fe 0cc8cc3b c9f61eab c4688f11 e9f64a2e"
#define P192_S_MINUS_T_Y "093e31d0 0fb78269 732b1bd2 a73c23cd d31745d0 523d816b"
#define P192_2S_X "30c5bc6b 8c7da253 54b373dc 14dd8a0e ba42d25a 3f6e6962"
#define P192_2S_Y "0dde14bc 4249a721 c407aedb f011e2dd bbcb2968 c9d889cf"
#define P192_SCALAR_D "a78a236d 60baec0c 5dd41b33 a542463a 8255391a f64c74ee"
#define P192_D_TIMES_S_X "1faee420 5a4f669d 2d0a8f25 e3bcec9a 62a69529 65bf6d31"
#define P192_D_TIMES_S_Y "5ff2cdfa 508a2581 89236708 7c696f17 9e7a4d7e 8260fb06"

#define P224_POINT_S_X "6eca814b a59a9308 43dc814e dd6c97da 95518df3 c6fdf16e 9a10bb5b"
#define P224_POINT_S_Y "ef4b497f 0963bc8b 6aec0ca0 f259b89c d8099414 7e05dc6b 64d7bf22"
#define P224_POINT_T_X "b72b25ae a5cb03fb 88d7e842 00296964 8e6ef23c 5d39ac90 3826bd6d"
#define P224_POINT_T_Y "c42a8a4d 34984f0b 71b5b409 1af7dceb 33ea729c 1a2dc8b4 34f10c34"
#define P224_S_PLUS_T_X "236f26d9 e84c2f7d 776b107b d478ee0a 6d2bcfca a2162afa e8d2fd15"
#define P224_S_PLUS_T_Y "e53cc0a7 904ce6c3 746f6a97 471297a0 b7d5cdf8 d536ae25 bb0fda70"
#define P224_S_MINUS_T_X "db4112bc c8f34d4f 0b36047b ca1054f3 61541385 2a793133 5210b332"
#define P224_S_MINUS_T_Y "90c6e830 4da48138 78c1540b 2396f411 facf787a 520a0ffb 55a8d961"
#define P224_2S_X "a9c96f21 17dee0f2 7ca56850 ebb46efa d8ee2685 2f165e29 cb5cdfc7"
#define P224_2S_Y "adf18c84 cf77ced4 d76d4930 417d9579 207840bf 49bfbf58 37dfdd7d"
#define P224_SCALAR_D "a78ccc30 eaca0fcc 8e36b2dd 6fbb03df 06d37f52 711e6363 aaf1d73b"
#define P224_D_TIMES_S_X "96a7625e 92a8d72b ff1113ab db95777e 736a14c6 fdaacc39 2702bca4"
#define P224_D_TIMES_S_Y "0f8e5702 942a3c5e 13cd2fd5 80191525 8b43dfad c70d15db ada3ed10"

#define P256_POINT_S_X "de2444be bc8d36e6 82edd27e 0f271508 617519b3 221a8fa0 b77cab39 89da97c9"
#define P256_POINT_S_Y "c093ae7f f36e5380 fc01a5aa d1e66659 702de80f 53cec576 b6350b24 3042a256"
#define P256_POINT_T_X "55a8b00f 8da1d44e 62f6b3b2 5316212e 39540dc8 61c89575 bb8cf92e 35e0986b"
#define P256_POINT_T_Y "5421c320 9c2d6c70 4835d82a c4c3dd90 f61a8a52 598b9e7a b656e9d8 c8b24316"
#define P256_S_PLUS_T_X "72b13dd4 354b6b81 745195e9 8cc5ba69 70349191 ac476bd4 553cf35a 545a067e"
#define P256_S_PLUS_T_Y "8d585cbb 2e1327d7 5241a8a1 22d7620d c33b1331 5aa5c9d4 6d013011 744ac264"
#define P256_S_MINUS_T_X "c09ce680 b251bb1d 2aad1dbf 6129deab 837419f8 f1c73ea1 3e7dc64a d6be6021"
#define P256_S_MINUS_T_Y "1a815bf7 00bd8833 6b2f9bad 4edab172 3414a022 fdf6c3f4 ce30675f b1975ef3"
#define P256_2S_X "7669e690 1606ee3b a1a8eef1 e0024c33 df6c22f3 b17481b8 2a860ffc db6127b0"
#define P256_2S_Y "fa878162 187a54f6 c39f6ee0 072f33de 389ef3ee cd03023d e10ca2c1 db61d0c7"
#define P256_SCALAR_D "c51e4753 afdec1e6 b6c6a5b9 92f43f8d d0c7a893 3072708b 6522468b 2ffb06fd"
#define P256_D_TIMES_S_X "51d08d5f 2d427888 2946d88d 83c97d11 e62becc3 cfc18bed acc89ba3 4eeca03f"
#define P256_D_TIMES_S_Y "75ee68eb 8bf626aa 5b673ab5 1f6e744e 06f8fcf8 a6c0cf30 35beca95 6a7b41d5"

#define P384_POINT_S_X "fba203b8 1bbd23f2 b3be971c c23997e1 ae4d89e6 9cb6f923 85dda827 68ada415 ebab4167 459da98e 62b1332d 1e73cb0e"
#define P384_POINT_S_Y "5ffedbae fdeba603 e7923e06 cdb5d0c6 5b223014 29293376 d5c6944e 3fa6259f 162b4788 de6987fd 59aed5e4 b5285e45"
#define P384_POINT_T_X "aacc0520 2e7fda6f c73d82f0 a6622052 7da8117e e8f8330e ad7d20ee 6f255f58 2d8bd38c 5a7f2b40 bcdb68ba 13d81051"
#define P384_POINT_T_Y "84009a26 3fefba7c 2c57cffa 5db3634d 286131af c0fca8d2 5afa22a7 b5dce0d9 470da892 33cee178 592f49b6 fecb5092"
#define P384_S_PLUS_T_X "12dc5ce7 acdfc584 4d939f40 b4df012e 68f865b8 9c3213ba 97090a24 7a2fc009 075cf471 cd2e85c4 89979b65 ee0b5eed"
#define P384_S_PLUS_T_Y "167312e5 8fe0c0af a248f285 4e3cddcb 557f983b 3189b67f 21eee013 41e7e9fe 67f6ee81 b36988ef a406945c 8804a4b0"
#define P384_S_MINUS_T_X "6afdaf8d a8b11c98 4cf177e5 51cee542 cda4ac2f 25cd522d 0cd710f8 8059c656 5aef78f6 b5ed6cc0 5a6666de f2a2fb59"
#define P384_S_MINUS_T_Y "7bed0e15 8ae8cc70 e847a603 47ca1548 c348decc 6309f48b 59bd5afc 9a9b804e 7f787617 8cb5a7eb 4f6940a9 c73e8e5e"
#define P384_2S_X "2a2111b1 e0aa8b2f c5a19755 16bc4d58 017ff96b 25e1bdff 3c229d5f ac3bacc3 19dcbec2 9f9478f4 2dee597b 4641504c"
#define P384_2S_Y "fa2e3d9d c84db895 4ce8085e f28d7184 fddfd134 4b4d4797 343af9b5 f9d83752 0b450f72 6443e411 4bd4e5bd b2f65ddd"
#define P384_SCALAR_D "a4ebcae5 a6659834 93ab3e62 6085a24c 104311a7 61b5a8fd ac052ed1 f111a5c4 4f76f456 59d2d111 a61b5fdd 97583480"
#define P384_D_TIMES_S_X "e4f77e7f feb7f095 8910e3a6 80d677a4 77191df1 66160ff7 ef6bb526 1f791aa7 b45e3e65 3d151b95 dad3d93c a0290ef2"
#define P384_D_TIMES_S_Y "ac7dee41 d8c5f4a7 d5836960 a773cfc1 376289d3 373f8cf7 417b0c62 07ac32e9 13856612 fc9ff2e3 57eb2ee0 5cf9667f"

#define P521_POINT_S_X "000001d5 c693f66c 08ed03ad 0f031f93 7443458f 601fd098 d3d0227b 4bf62873 af50740b 0bb84aa1 57fc847b cf8dc16a 8b2b8bfd 8e2d0a7d 39af04b0 89930ef6 dad5c1b4"
#define P521_POINT_S_Y "00000144 b7770963 c63a3924 8865ff36 b074151e ac33549b 224af5c8 664c5401 2b818ed0 37b2b7c1 a63ac89e baa11e07 db89fcee 5b556e49 764ee3fa 66ea7ae6 1ac01823"
#define P521_POINT_T_X "000000f4 11f2ac2e b971a267 b80297ba 67c322db a4bb21ce c8b70073 bf88fc1c a5fde3ba 09e5df6d 39acb2c0 762c03d7 bc224a3e 197feaf7 60d63240 06fe3be9 a548c7d5"
#define P521_POINT_T_Y "000001fd f842769c 707c93c6 30df6d02 eff399a0 6f1b36fb 9684f0b3 73ed0648 89629abb 92b1ae32 8fdb4553 42683849 43f0e922 2afe0325 9b32274d 35d1b958 4c65e305"
#define P521_S_PLUS_T_X "00000126 4ae115ba 9cbc2ee5 6e6f0059 e24b52c8 04632160 2c59a339 cfb757c8 9a59c358 a9a8e1f8 6d384b3f 3b255ea3 f73670c6 dc9f45d4 6b6a196d c37bbe0f 6b2dd9e9"
#define P521_S_PLUS_T_Y "00000062 a9c72b8f 9f88a271 690bfa01 7a6466c3 1b9cadc2 fc544744 aeb81707 2349cfdd c5ad0e81 b03f1897 bd9c8c6e fbdf6823 7dc3bb00 445979fb 373b20c9 a967ac55"
#define P521_S_MINUS_T_X "00000129 2cb58b17 95ba4770 63fef7cd 22e42c20 f57ae94c eaad86e0 d21ff229 18b0dd3b 076d63be 253de24b c20c6da290fa54d8 3771a225 deecf914 9f79a8e6 14c3c4cd"
#define P521_S_MINUS_T_Y "00000169 5e3821e7 2c7cacaa dcf62909 cd83463a 21c6d033 93c527c6 43b36239 c46af117 ab7c7ad1 9a4c8cf0 ae95ed51 72988546 1aa2ce27 00a6365b ca3733d2 920b2267"
#define P521_2S_X "00000128 79442f24 50c119e7 119a5f73 8be1f1eb a9e9d7c6 cf41b325 d9ce6d64 3106e9d6 1124a91a 96bcf201 305a9dee 55fa7913 6dc70083 1e54c3ca 4ff2646b d3c36bc6"
#define P521_2S_Y "00000198 64a8b885 5c2479cb efe375ae 553e2393 271ed36f adfc4494 fc0583f6 bd035988 96f39854 abeae5f9 a6515a02 1e2c0eef 139e71de 610143f5 3382f410 4dccb543"
#define P521_SCALAR_D "000001eb 7f81785c 9629f136 a7e8f8c6 74957109 73555411 1a2a866f a5a16669 9419bfa9 936c78b6 2653964d f0d6da94 0a695c72 94d41b2d 6600de6d fcf0edcf c89fdcb1"
#define P521_D_TIMES_S_X "00000091 b15d09d0 ca0353f8 f96b93cd b13497b0 a4bb582a e9ebefa3 5eee61bf 7b7d041b 8ec34c6c 00c0c067 1c4ae063 318fb75b e87af4fe 859608c9 5f0ab477 4f8c95bb"
#define P521_D_TIMES_S_Y "00000130 f8f8b5e1 abb4dd94 f6baaf65 4a2d5810 411e77b7 423965e0 c7fd79ec 1ae563c2 07bd255e e9828eb7 a03fed56 5240d2cc 80ddd2ce cbb2eb50 f0951f75 ad87977f"


class CurvesNISTTest : public CppUnit::TestFixture {
	CPPUNIT_TEST_SUITE(CurvesNISTTest);

	CPPUNIT_TEST(p192Addition);
	CPPUNIT_TEST(p192Subtraction);
	CPPUNIT_TEST(p192Doubling);
	CPPUNIT_TEST(p192Multiplication);
	// Tests if p^ord=e
	CPPUNIT_TEST(p192Order);

	CPPUNIT_TEST(p224Addition);
	CPPUNIT_TEST(p224Subtraction);
	CPPUNIT_TEST(p224Doubling);
	CPPUNIT_TEST(p224Multiplication);
	CPPUNIT_TEST(p224Order);

	CPPUNIT_TEST(p256Addition);
	CPPUNIT_TEST(p256Subtraction);
	CPPUNIT_TEST(p256Doubling);
	CPPUNIT_TEST(p256Multiplication);
	CPPUNIT_TEST(p256Order);

	CPPUNIT_TEST(p384Addition);
	CPPUNIT_TEST(p384Subtraction);
	CPPUNIT_TEST(p384Doubling);
	CPPUNIT_TEST(p384Multiplication);
	CPPUNIT_TEST(p384Order);

	CPPUNIT_TEST(p521Addition);
	CPPUNIT_TEST(p521Subtraction);
	CPPUNIT_TEST(p521Doubling);
	CPPUNIT_TEST(p521Multiplication);
	CPPUNIT_TEST(p521Order);

	CPPUNIT_TEST_SUITE_END();

public:
	CurvesNISTTest();

	void setUp();
	void tearDown();

	void p192Addition();
	void p192Subtraction();
	void p192Doubling();
	void p192Multiplication();
	void p192Order();

	void p224Addition();
	void p224Subtraction();
	void p224Doubling();
	void p224Multiplication();
	void p224Order();

	void p256Addition();
	void p256Subtraction();
	void p256Doubling();
	void p256Multiplication();
	void p256Order();

	void p384Addition();
	void p384Subtraction();
	void p384Doubling();
	void p384Multiplication();
	void p384Order();

	void p521Addition();
	void p521Subtraction();
	void p521Doubling();
	void p521Multiplication();
	void p521Order();

private:
	ECPrime *curveP192, *curveP224, *curveP256, *curveP384, *curveP521;
	Coordinate p192S, p192T, p224S, p224T, p256S, p256T, p384S, p384T, p521S, p521T;
	mpz_class p192d, p224d, p256d, p384d, p521d;

};

#endif /* CURVESNISTTEST_H_ */
