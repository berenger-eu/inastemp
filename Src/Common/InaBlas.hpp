///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECINTERFACE_HPP
#define INAVECINTERFACE_HPP

#include <cstdio>
#include <vector>

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


    inline void setScalar(RealType* ptr, const RealType inVal, const long int nbValues){
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

    inline void VecMultScalar(RealType* ptr, const RealType scalar, const long int nbValues){
        VecType vec_scal = scalar;
        VecType vec_ptr;
        VecType res;

        const long int nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr) & (Alignement-1)) == 0) {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr.setFromAlignedArray(&ptr[idx])*vec_scal;
                res.storeInAlignedArray(&ptr[idx]);
                res=vec_ptr.setFromAlignedArray(&ptr[idx+VecType::GetVecLength()])*vec_scal;
                res.storeInAlignedArray(&ptr[idx+VecType::GetVecLength()]);
                res=vec_ptr.setFromAlignedArray(&ptr[idx+2*VecType::GetVecLength()])*vec_scal;
                res.storeInAlignedArray(&ptr[idx+2*VecType::GetVecLength()]);
                res=vec_ptr.setFromAlignedArray(&ptr[idx+3*VecType::GetVecLength()])*vec_scal;
                res.storeInAlignedArray(&ptr[idx+3*VecType::GetVecLength()]);
            }
        }
        else {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr.setFromArray(&ptr[idx])*vec_scal;
                res.storeInArray(&ptr[idx]);
                res=vec_ptr.setFromArray(&ptr[idx+VecType::GetVecLength()])*vec_scal;
                res.storeInArray(&ptr[idx+VecType::GetVecLength()]);
                res=vec_ptr.setFromArray(&ptr[idx+2*VecType::GetVecLength()])*vec_scal;
                res.storeInArray(&ptr[idx+2*VecType::GetVecLength()]);
                res=vec_ptr.setFromArray(&ptr[idx+3*VecType::GetVecLength()])*vec_scal;
                res.storeInArray(&ptr[idx+3*VecType::GetVecLength()]);
            }
        }

        for(long int idx = nbValuesRounded4 ; idx < nbValues ; ++idx){
            ptr[idx] *= scalar;
        }

    }


    inline void VecAddScalar(RealType* ptr, const RealType scalar, const long int nbValues){
        VecType vec_scal = scalar;
        VecType vec_ptr;
        VecType res;

        const long int nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr) & (Alignement-1)) == 0) {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr.setFromAlignedArray(&ptr[idx])+vec_scal;
                res.storeInAlignedArray(&ptr[idx]);
                res=vec_ptr.setFromAlignedArray(&ptr[idx+VecType::GetVecLength()])+vec_scal;
                res.storeInAlignedArray(&ptr[idx+VecType::GetVecLength()]);
                res=vec_ptr.setFromAlignedArray(&ptr[idx+2*VecType::GetVecLength()])+vec_scal;
                res.storeInAlignedArray(&ptr[idx+2*VecType::GetVecLength()]);
                res=vec_ptr.setFromAlignedArray(&ptr[idx+3*VecType::GetVecLength()])+vec_scal;
                res.storeInAlignedArray(&ptr[idx+3*VecType::GetVecLength()]);
            }
        }
        else {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr.setFromArray(&ptr[idx])+vec_scal;
                res.storeInArray(&ptr[idx]);
                res=vec_ptr.setFromArray(&ptr[idx+VecType::GetVecLength()])+vec_scal;
                res.storeInArray(&ptr[idx+VecType::GetVecLength()]);
                res=vec_ptr.setFromArray(&ptr[idx+2*VecType::GetVecLength()])+vec_scal;
                res.storeInArray(&ptr[idx+2*VecType::GetVecLength()]);
                res=vec_ptr.setFromArray(&ptr[idx+3*VecType::GetVecLength()])+vec_scal;
                res.storeInArray(&ptr[idx+3*VecType::GetVecLength()]);
            }
        }

        for(long int idx = nbValuesRounded4 ; idx < nbValues ; ++idx){
            ptr[idx] += scalar;
        }

    }


    inline RealType* MultTermToTerm(RealType* ptr1, RealType* ptr2, const long int nbValues){
        VecType vec_ptr1, vec_ptr2, res;

        alignas(512) char buff2[VecType::GetVecLength()*sizeof(RealType)+1];
        RealType* ptr_res = reinterpret_cast<RealType*>(&buff2[1]);

        const long int nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr1) & (Alignement-1)) == 0) {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx])*vec_ptr2.setFromAlignedArray(&ptr2[idx]);
                res.storeInAlignedArray(&ptr_res[idx]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+VecType::GetVecLength()])*vec_ptr2.setFromAlignedArray(&ptr2[idx+VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+VecType::GetVecLength()]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+2*VecType::GetVecLength()])*vec_ptr2.setFromAlignedArray(&ptr2[idx+3*VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+2*VecType::GetVecLength()]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+3*VecType::GetVecLength()])*vec_ptr2.setFromAlignedArray(&ptr2[idx+4*VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+3*VecType::GetVecLength()]);
            }
        }
        else {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromArray(&ptr1[idx])*vec_ptr2.setFromArray(&ptr2[idx]);
                res.storeInArray(&ptr_res[idx]);
                res=vec_ptr1.setFromArray(&ptr1[idx+VecType::GetVecLength()])*vec_ptr2.setFromArray(&ptr2[idx+VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+VecType::GetVecLength()]);
                res=vec_ptr1.setFromArray(&ptr1[idx+2*VecType::GetVecLength()])*vec_ptr2.setFromArray(&ptr2[idx+3*VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+2*VecType::GetVecLength()]);
                res=vec_ptr1.setFromArray(&ptr1[idx+3*VecType::GetVecLength()])*vec_ptr2.setFromArray(&ptr2[idx+4*VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+3*VecType::GetVecLength()]);
            }
        }

        for(long int idx = nbValuesRounded4 ; idx < nbValues ; ++idx){
            ptr_res[idx] = ptr1[idx] * ptr2[idx];
        }

        return ptr_res;

    }


    inline RealType* AddTermToTerm(RealType* ptr1, RealType* ptr2, const long int nbValues){
        VecType vec_ptr1, vec_ptr2, res;

        alignas(512) char buff2[VecType::GetVecLength()*sizeof(RealType)+1];
        RealType* ptr_res = reinterpret_cast<RealType*>(&buff2[1]);

        const long int nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr1) & (Alignement-1)) == 0) {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx])+vec_ptr2.setFromAlignedArray(&ptr2[idx]);
                res.storeInAlignedArray(&ptr_res[idx]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+VecType::GetVecLength()])+vec_ptr2.setFromAlignedArray(&ptr2[idx+VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+VecType::GetVecLength()]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+2*VecType::GetVecLength()])+vec_ptr2.setFromAlignedArray(&ptr2[idx+3*VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+2*VecType::GetVecLength()]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+3*VecType::GetVecLength()])+vec_ptr2.setFromAlignedArray(&ptr2[idx+4*VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+3*VecType::GetVecLength()]);
            }
        }
        else {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromArray(&ptr1[idx])+vec_ptr2.setFromArray(&ptr2[idx]);
                res.storeInArray(&ptr_res[idx]);
                res=vec_ptr1.setFromArray(&ptr1[idx+VecType::GetVecLength()])+vec_ptr2.setFromArray(&ptr2[idx+VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+VecType::GetVecLength()]);
                res=vec_ptr1.setFromArray(&ptr1[idx+2*VecType::GetVecLength()])+vec_ptr2.setFromArray(&ptr2[idx+3*VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+2*VecType::GetVecLength()]);
                res=vec_ptr1.setFromArray(&ptr1[idx+3*VecType::GetVecLength()])+vec_ptr2.setFromArray(&ptr2[idx+4*VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+3*VecType::GetVecLength()]);
            }
        }

        for(long int idx = nbValuesRounded4 ; idx < nbValues ; ++idx){
            ptr_res[idx] = ptr1[idx] + ptr2[idx];
        }

        return ptr_res;

    }

};




#endif // INABLAS_HPP
