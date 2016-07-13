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
inline typename VecType::VecRawType FastPow(typename VecType::VecRawType base, size_t power){
    typename VecType::VecRawType res = VecType(1.).getVec();

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
