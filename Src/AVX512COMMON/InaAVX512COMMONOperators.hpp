///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAAVX512COMMONOPERATORS_HPP
#define INAAVX512COMMONOPERATORS_HPP

#include "InastempConfig.h"

#ifndef INASTEMP_USE_AVX512COMMON
#error InaAVX512COMMONOperators is included but AVX512COMMON is not enable in the configuration
#endif

#include <immintrin.h>
#include <cmath>

#ifdef INASTEMP_USE_AVX512COMMON_OPERATORS

//Side effect operators DOUBLE
inline __m512d& operator+=(__m512d& a, const __m512d& b) {
    return (a = _mm512_add_pd(a, b));
}

inline __m512d& operator-=(__m512d& a, const __m512d& b) {
    return (a = _mm512_sub_pd(a, b));
}

inline __m512d& operator*=(__m512d& a, const __m512d& b) {
    return (a = _mm512_mul_pd(a, b));
}

inline __m512d& operator/=(__m512d& a, const __m512d& b) {
    return (a = _mm512_div_pd(a, b));
}

//No side effect operators DOUBLE
inline __m512d operator+(const __m512d& a, const __m512d& b) {
    return _mm512_add_pd(a, b);
}

inline __m512d operator-(const __m512d& a, const __m512d& b) {
    return _mm512_sub_pd(a, b);
}

inline __m512d operator*(const __m512d& v1, const __m512d& v2) {
    return _mm512_mul_pd(v1, v2);
}

inline __m512d operator/(const __m512d& v1, const __m512d& v2) {
    return _mm512_div_pd(v1, v2);
}

//Side effect operators SINGLE
inline __m512& operator+=(__m512& a, const __m512& b) {
    return (a = _mm512_add_ps(a, b));
}

inline __m512& operator-=(__m512& a, const __m512& b) {
    return (a = _mm512_sub_ps(a, b));
}

inline __m512& operator*=(__m512& a, const __m512& b) {
    return (a = _mm512_mul_ps(a, b));
}

inline __m512& operator/=(__m512& a, const __m512& b) {
    return (a = _mm512_div_ps(a, b));
}

//No side effect operators SINGLE
inline __m512 operator+(const __m512& a, const __m512& b) {
    return _mm512_add_ps(a, b);
}

inline __m512 operator-(const __m512& a, const __m512& b) {
    return _mm512_sub_ps(a, b);
}

inline __m512 operator*(const __m512& v1, const __m512& v2) {
    return _mm512_mul_ps(v1, v2);
}

inline __m512 operator/(const __m512& v1, const __m512& v2) {
    return _mm512_div_ps(v1, v2);
}

#endif

#endif
