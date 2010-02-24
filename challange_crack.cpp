#include "challange_crack.h"
#include "hcp.h"
#include "elgamal.h"

#define MERSENNE_PRIME_NUMBER 13
ECPrime candidate_elliptic_curve(int mersenne_prime_number, int D, string ECPoint){
    ECPrime curve = ECPrime::normalizedCurveFromDiscriminantAndPrime(D, mersenne_prime(mersenne_prime_number));
    curve.set_point_compressed(ECPoint);
    return curve;
}

void try_all_curves_on_cipher(string ECPoint, string d_string, string cipher){
    #define LOWER_LIMIT -400
    HCP temp;
    mpz_class d;
    mpz_set_str(d.get_mpz_t(), d_string.c_str(), 10);
    for (int D=-1; D>LOWER_LIMIT; D--){
        if ((temp.H.find(D)) != temp.H.end()){//this gives something to work with
            ECPrime curve = candidate_elliptic_curve(MERSENNE_PRIME_NUMBER, D, ECPoint);
            ECC_ElGamal elgamal(curve);
            elgamal.set_private_key(d);
            cout << elgamal.decrypt(cipher) << endl;
        }
    }
}
