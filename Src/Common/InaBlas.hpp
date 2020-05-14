///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECINTERFACE_HPP
#define INAVECINTERFACE_HPP

#include <cstdio>

template < class VecType >
class InaBlas : public VecType {
public:
    using RealType = typename VecType::RealType;
    static const int Alignement= VecType::Alignement;

    inline InaBlas(){}


    inline void setZero(RealType* ptr, const long int nbValues){
        VecType zero = VecType::GetZero();
        const long int nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr) & (Alignement-1)) == 0) {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                 zero.storeInAlignedArray(&ptr[idx]);
                 zero.storeInAlignedArray(&ptr[idx+VecType::GetVecLength()]);
                 zero.storeInAlignedArray(&ptr[idx+2*VecType::GetVecLength()]);
                 zero.storeInAlignedArray(&ptr[idx+3*VecType::GetVecLength()]);
             }
        }
        else {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                 zero.storeInArray(&ptr[idx]);
                 zero.storeInArray(&ptr[idx+VecType::GetVecLength()]);
                 zero.storeInArray(&ptr[idx+2*VecType::GetVecLength()]);
                 zero.storeInArray(&ptr[idx+3*VecType::GetVecLength()]);
             }
        }

         for(long int idx = nbValuesRounded4 ; idx < nbValues ; ++idx){
             ptr[idx] = 0;
         }

    }


    inline void setScalar(RealType* ptr, const long int nbValues, const RealType inVal){
        VecType vec = inVal;

        const long int nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr) & (Alignement-1)) == 0) {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                 vec.storeInAlignedArray(&ptr[idx]);
                 vec.storeInAlignedArray(&ptr[idx+VecType::GetVecLength()]);
                 vec.storeInAlignedArray(&ptr[idx+2*VecType::GetVecLength()]);
                 vec.storeInAlignedArray(&ptr[idx+3*VecType::GetVecLength()]);
             }
        }
        else {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                 vec.storeInArray(&ptr[idx]);
                 vec.storeInArray(&ptr[idx+VecType::GetVecLength()]);
                 vec.storeInArray(&ptr[idx+2*VecType::GetVecLength()]);
                 vec.storeInArray(&ptr[idx+3*VecType::GetVecLength()]);
             }
        }

         for(long int idx = nbValuesRounded4 ; idx < nbValues ; ++idx){
             ptr[idx] = inVal;
         }

    }

 };




#endif // INABLAS_HPP
