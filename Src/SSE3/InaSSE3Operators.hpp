///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INASSE3OPERATORS_HPP
#define INASSE3OPERATORS_HPP

#include "InastempConfig.h"

#ifndef INASTEMP_USE_SSE3
#error InaSSE3Operators is included but SSE3 is not enable in the configuration
#endif

#include <emmintrin.h>
#include <cmath>

#ifdef INASTEMP_USE_SSE3_OPERATORS

inline __m128d& operator+=(__m128d& v1, const __m128d& v2) {
    return (v1 = _mm_add_pd(v1, v2));
}

inline __m128d& operator-=(__m128d& v1, const __m128d& v2) {
    return (v1 = _mm_sub_pd(v1, v2));
}

inline __m128d& operator*=(__m128d& v1, const __m128d& v2) {
    return (v1 = _mm_mul_pd(v1, v2));
}

inline __m128d& operator/=(__m128d& v1, const __m128d& v2) {
    return (v1 = _mm_div_pd(v1, v2));
}

inline __m128d operator+(const __m128d& v1, const __m128d& v2) {
    return _mm_add_pd(v1, v2);
}

inline __m128d operator-(const __m128d& v1, const __m128d& v2) {
    return _mm_sub_pd(v1, v2);
}

inline __m128d operator*(const __m128d& v1, const __m128d& v2) {
    return _mm_mul_pd(v1, v2);
}

inline __m128d operator/(const __m128d& v1, const __m128d& v2) {
    return _mm_div_pd(v1, v2);
}

inline __m128& operator+=(__m128& v1, const __m128& v2) {
    return (v1 = _mm_add_ps(v1, v2));
}

inline __m128& operator-=(__m128& v1, const __m128& v2) {
    return (v1 = _mm_sub_ps(v1, v2));
}

inline __m128& operator*=(__m128& v1, const __m128& v2) {
    return (v1 = _mm_mul_ps(v1, v2));
}

inline __m128& operator/=(__m128& v1, const __m128& v2) {
    return (v1 = _mm_div_ps(v1, v2));
}

inline __m128 operator+(const __m128& v1, const __m128& v2) {
    return _mm_add_ps(v1, v2);
}

inline __m128 operator-(const __m128& v1, const __m128& v2) {
    return _mm_sub_ps(v1, v2);
}

inline __m128 operator*(const __m128& v1, const __m128& v2) {
    return _mm_mul_ps(v1, v2);
}

inline __m128 operator/(const __m128& v1, const __m128& v2) {
    return _mm_div_ps(v1, v2);
}

#endif

#endif
