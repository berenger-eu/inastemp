///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef UTILS_HPP
#define UTILS_HPP

// for std::size_t
#include <cstddef>
#include <cstring>
#include <type_traits>

namespace InaUtils {

template < class TypeSrc, class TypeDest >
inline TypeDest ConvertBits(TypeSrc src) {
    static_assert(sizeof(TypeSrc) == sizeof(TypeDest), "Bit conversion requires equal sizes");
    static_assert(std::is_trivially_copyable<TypeSrc>::value &&
                  std::is_trivially_copyable<TypeDest>::value,
                  "Bit conversion requires trivially copyable types");
    TypeDest dest;
    std::memcpy(&dest, &src, sizeof(dest));
    return dest;
}

template <class VecType>
inline VecType FastPow(VecType base, std::size_t power){
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

inline std::size_t FastPowNbMul(std::size_t power){
    std::size_t nbMul = 0;

    while(power){
        if(1 & power){
            nbMul += 1;
        }
        nbMul += 1;
        power >>= 1;
    }

    return nbMul;
}

}

#endif // UTILS_HPP
