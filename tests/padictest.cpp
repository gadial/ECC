/*
 * padictest.cpp
 *
 *  Created on: Jan 23, 2010
 *      Author: bhess
 */

#include "padictest.h"

void Padictest::setUp() {
	M = new Poly(8);
	Poly M(8);
	M.set_coeff(0, 1);
	M.set_coeff(2, 1);
	M.set_coeff(3, 1);
	M.set_coeff(4, 1);
	M.set_coeff(8, 1);

	ad = new Adicops(M);
	ad->set_teichmuller_modulus(M, 10);


}

void Padictest::tearDown() {
	delete M;
	delete ad;
}

void Padictest::teichmueller_mod() {
	Poly mod = ad->get_mod();
	Poly v(8);
	v.set_coeff(8, 1);
	v.set_coeff(7, 644);
	v.set_coeff(6, 842);
	v.set_coeff(5, 134);
	v.set_coeff(4, 523);
	v.set_coeff(3, 21);
	v.set_coeff(2, 1019);
	v.set_coeff(1, 562);
	v.set_coeff(0, 1);
	CPPUNIT_ASSERT(v == mod);
}

void Padictest::fast_division_with_remainder() {
	Poly A(14);
	A.set_coeff(14, 559);
	A.set_coeff(13, 781);
	A.set_coeff(12, 763);
	A.set_coeff(11, 684);
	A.set_coeff(10, 133);
	A.set_coeff(9, 375);
	A.set_coeff(8, 922);
	A.set_coeff(7, 776);
	A.set_coeff(6, 452);
	A.set_coeff(5, 214);
	A.set_coeff(4, 313);
	A.set_coeff(3, 148);
	A.set_coeff(2, 646);
	A.set_coeff(1, 428);
	A.set_coeff(0, 168);

	Poly q = ad->poly_division_rem(A, ad->get_mod(), 10);
	Poly r = ad->poly_remainder(A, 10);

	Poly refQ(6);
	refQ.set_coeff(6, 559);
	refQ.set_coeff(5, 209);
	refQ.set_coeff(4, 673);
	refQ.set_coeff(3, 420);
	refQ.set_coeff(2, 768);
	refQ.set_coeff(1, 755);
	refQ.set_coeff(0, 337);
	CPPUNIT_ASSERT(refQ == q);

	Poly refR(7);
	refR.set_coeff(7, 428);
	refR.set_coeff(6, 728);
	refR.set_coeff(5, 240);
	refR.set_coeff(4, 294);
	refR.set_coeff(3, 10);
	refR.set_coeff(2, 165);
	refR.set_coeff(1, 743);
	refR.set_coeff(0, 855);
	CPPUNIT_ASSERT(refR == r);
}

void Padictest::inverse() {
	Poly A = Poly(7);
	A.set_coeff(7, 982);
	A.set_coeff(6, 303);
	A.set_coeff(5, 724);
	A.set_coeff(4, 458);
	A.set_coeff(3, 918);
	A.set_coeff(2, 423);
	A.set_coeff(1, 650);
	A.set_coeff(0, 591);
	//A.print();

	Poly invA = ad->get_inverse(A, 10);

	Poly veriInvA(7);
	veriInvA.set_coeff(7, 854);
	veriInvA.set_coeff(6, 373);
	veriInvA.set_coeff(5, 760);
	veriInvA.set_coeff(4, 132);
	veriInvA.set_coeff(3, 863);
	veriInvA.set_coeff(2, 697);
	veriInvA.set_coeff(1, 321);
	veriInvA.set_coeff(0, 60);

	CPPUNIT_ASSERT(invA == veriInvA);
}

void Padictest::invsqrt() {
	Poly AA(7);
	AA.set_coeff(7, 823);
	AA.set_coeff(6, 707);
	AA.set_coeff(5, 860);
	AA.set_coeff(4, 387);
	AA.set_coeff(3, 663);
	AA.set_coeff(2, 183);
	AA.set_coeff(1, 12);
	AA.set_coeff(0, 354);
	//AA.print();

	Poly AAaprox(7);
	AAaprox.set_coeff(7, 2);
	AAaprox.set_coeff(6, 1);
	AAaprox.set_coeff(3, 3);
	AAaprox.set_coeff(2, 1);
	AAaprox.set_coeff(1, 1);

	Poly insq = ad->get_invsqrt(AA, AAaprox, 9);

	Poly veriInsq(7);
	veriInsq.set_coeff(7, 342);
	veriInsq.set_coeff(6, 373);
	veriInsq.set_coeff(5, 248);
	veriInsq.set_coeff(4, 132);
	veriInsq.set_coeff(3, 351);
	veriInsq.set_coeff(2, 185);
	veriInsq.set_coeff(1, 321);
	veriInsq.set_coeff(0, 60);

	CPPUNIT_ASSERT(veriInsq == insq);
}
