///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef UTILS_HPP
#define UTILS_HPP

namespace InaUtils {

template < class TypeSrc, class TypeDest >
inline TypeDest ConvertBits(TypeSrc src) {
    union InaBits {
        TypeSrc src;
        TypeDest dest;
    };
    InaBits bits;
    bits.src = src;
    return bits.dest;
}

template <class VecType>
inline VecType FastPow(VecType base, size_t power){
    VecType res = VecType(1.);

    while(power){
        if(1 & power){
            res *= base;
        }
        base *= base;
        power >>= 1;
    }

    return res;
}

}

#endif // UTILS_HPP
