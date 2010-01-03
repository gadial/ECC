#include "zp_int.h"

zp_int& zp_int::operator+=(const zp_int& rhs){
    val += rhs.val;
    return normalize();
}

zp_int& zp_int::operator-=(const zp_int& rhs){
    val -= rhs.val;
    return normalize();
}

zp_int& zp_int::operator*=(const zp_int& rhs){
    val *= rhs.val;
    return normalize();
}

zp_int& zp_int::operator/=(const zp_int& rhs){
    val *= rhs.inverse().val;
    return normalize();
}

zp_int& zp_int::invert(){
    if (val == 0)
        throw "Division by 0";
    mpz_invert(val.get_mpz_t(), val.get_mpz_t(),p.get_mpz_t());
    return normalize();
}

zp_int zp_int::inverse() const{
    zp_int temp = *this;
    return temp.invert();
}

zp_int& zp_int::normalize(){
    if (p == 0)
        return *this;
    val %= p;
    if (val < 0)
        val += p;
    return *this;
}
bool zp_int::is_equal(const zp_int& rhs) const{
    return (val == rhs.val);
}

zp_int& zp_int::operator=(const zp_int& rhs){
    val = rhs.val;
    p = rhs.p;
}

zp_int operator+(const zp_int& lhs, const zp_int& rhs){
    zp_int temp = lhs;
    return (temp += rhs);
}
zp_int operator-(const zp_int& lhs, const zp_int& rhs){
    zp_int temp = lhs;
    return (temp -= rhs);
}
zp_int operator*(const zp_int& lhs, const zp_int& rhs){
    zp_int temp = lhs;
    return (temp *= rhs);
}
zp_int operator/(const zp_int& lhs, const zp_int& rhs){
    zp_int temp = lhs;
    return (temp /= rhs);
}

bool operator==(const zp_int& lhs, const zp_int& rhs){
    return (lhs.is_equal(rhs));
}

bool operator!=(const zp_int& lhs, const zp_int& rhs){
    return (!lhs.is_equal(rhs));
}