///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECINTERFACE_HPP
#define INAVECINTERFACE_HPP

#include <cstdio>
#include <cstring>
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
        VecType vect = inVal;

        const long int nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr) & (Alignement-1)) == 0) {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                vect.storeInAlignedArray(&ptr[idx]);
                vect.storeInAlignedArray(&ptr[idx+VecType::GetVecLength()]);
                vect.storeInAlignedArray(&ptr[idx+2*VecType::GetVecLength()]);
                vect.storeInAlignedArray(&ptr[idx+3*VecType::GetVecLength()]);
            }
        }
        else {
            for(long int idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                vect.storeInArray(&ptr[idx]);
                vect.storeInArray(&ptr[idx+VecType::GetVecLength()]);
                vect.storeInArray(&ptr[idx+2*VecType::GetVecLength()]);
                vect.storeInArray(&ptr[idx+3*VecType::GetVecLength()]);
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


    inline RealType* MultTermToTerm(RealType* ptr1, RealType* ptr2, const unsigned long nbValues){
        VecType vec_ptr1, vec_ptr2, res;

        RealType* ptr_res = new RealType[nbValues];

        const unsigned long nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr1) & (Alignement-1)) == 0) {
            for(unsigned long idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx])*vec_ptr2.setFromAlignedArray(&ptr2[idx]);
                res.storeInAlignedArray(&ptr_res[idx]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+VecType::GetVecLength()])*vec_ptr2.setFromAlignedArray(&ptr2[idx+VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+VecType::GetVecLength()]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+2*VecType::GetVecLength()])*vec_ptr2.setFromAlignedArray(&ptr2[idx+2*VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+2*VecType::GetVecLength()]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+3*VecType::GetVecLength()])*vec_ptr2.setFromAlignedArray(&ptr2[idx+3*VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+3*VecType::GetVecLength()]);
            }
        }
        else {
            for(unsigned long idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromArray(&ptr1[idx])*vec_ptr2.setFromArray(&ptr2[idx]);
                res.storeInArray(&ptr_res[idx]);
                res=vec_ptr1.setFromArray(&ptr1[idx+VecType::GetVecLength()])*vec_ptr2.setFromArray(&ptr2[idx+VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+VecType::GetVecLength()]);
                res=vec_ptr1.setFromArray(&ptr1[idx+2*VecType::GetVecLength()])*vec_ptr2.setFromArray(&ptr2[idx+2*VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+2*VecType::GetVecLength()]);
                res=vec_ptr1.setFromArray(&ptr1[idx+3*VecType::GetVecLength()])*vec_ptr2.setFromArray(&ptr2[idx+3*VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+3*VecType::GetVecLength()]);
            }
        }

        for(unsigned long idx = nbValuesRounded4 ; idx < nbValues ; ++idx){
            ptr_res[idx] = ptr1[idx] * ptr2[idx];
        }

        return ptr_res;

    }


    inline RealType* AddTermToTerm(RealType* ptr1, RealType* ptr2, const unsigned long nbValues){
        VecType vec_ptr1, vec_ptr2, res;

        RealType* ptr_res = new RealType[nbValues];

        const unsigned long nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr1) & (Alignement-1)) == 0) {
            for(unsigned long idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx])+vec_ptr2.setFromAlignedArray(&ptr2[idx]);
                res.storeInAlignedArray(&ptr_res[idx]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+VecType::GetVecLength()])+vec_ptr2.setFromAlignedArray(&ptr2[idx+VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+VecType::GetVecLength()]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+2*VecType::GetVecLength()])+vec_ptr2.setFromAlignedArray(&ptr2[idx+2*VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+2*VecType::GetVecLength()]);
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+3*VecType::GetVecLength()])+vec_ptr2.setFromAlignedArray(&ptr2[idx+3*VecType::GetVecLength()]);
                res.storeInAlignedArray(&ptr_res[idx+3*VecType::GetVecLength()]);
            }
        }
        else {
            for(unsigned long idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromArray(&ptr1[idx])+vec_ptr2.setFromArray(&ptr2[idx]);
                res.storeInArray(&ptr_res[idx]);
                res=vec_ptr1.setFromArray(&ptr1[idx+VecType::GetVecLength()])+vec_ptr2.setFromArray(&ptr2[idx+VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+VecType::GetVecLength()]);
                res=vec_ptr1.setFromArray(&ptr1[idx+2*VecType::GetVecLength()])+vec_ptr2.setFromArray(&ptr2[idx+2*VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+2*VecType::GetVecLength()]);
                res=vec_ptr1.setFromArray(&ptr1[idx+3*VecType::GetVecLength()])+vec_ptr2.setFromArray(&ptr2[idx+3*VecType::GetVecLength()]);
                res.storeInArray(&ptr_res[idx+3*VecType::GetVecLength()]);
            }
        }

        for(unsigned long idx = nbValuesRounded4 ; idx < nbValues ; ++idx){
            ptr_res[idx] = ptr1[idx] + ptr2[idx];
        }

        return ptr_res;

    }


    inline RealType ScalarProduct(RealType* ptr1, RealType* ptr2, const unsigned long nbValues){
        VecType vec_ptr1, vec_ptr2, res;
        RealType scalarProduct = RealType(0);

        const unsigned long nbValuesRounded4 = nbValues - (nbValues%(4*VecType::GetVecLength()));

        if((std::ptrdiff_t (ptr1) & (Alignement-1)) == 0) {
            for(unsigned long idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx])*vec_ptr2.setFromAlignedArray(&ptr2[idx]);
                scalarProduct+=res.horizontalSum();
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+VecType::GetVecLength()])*vec_ptr2.setFromAlignedArray(&ptr2[idx+VecType::GetVecLength()]);
                scalarProduct+=res.horizontalSum();
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+2*VecType::GetVecLength()])*vec_ptr2.setFromAlignedArray(&ptr2[idx+2*VecType::GetVecLength()]);
                scalarProduct+=res.horizontalSum();
                res=vec_ptr1.setFromAlignedArray(&ptr1[idx+3*VecType::GetVecLength()])*vec_ptr2.setFromAlignedArray(&ptr2[idx+3*VecType::GetVecLength()]);
                scalarProduct+=res.horizontalSum();
            }
        }
        else {
            for(unsigned long idx = 0 ; idx < nbValuesRounded4 ; idx += 4*VecType::GetVecLength()){
                res=vec_ptr1.setFromArray(&ptr1[idx])*vec_ptr2.setFromArray(&ptr2[idx]);
                scalarProduct+=res.horizontalSum();
                res=vec_ptr1.setFromArray(&ptr1[idx+VecType::GetVecLength()])*vec_ptr2.setFromArray(&ptr2[idx+VecType::GetVecLength()]);
                scalarProduct+=res.horizontalSum();
                res=vec_ptr1.setFromArray(&ptr1[idx+2*VecType::GetVecLength()])*vec_ptr2.setFromArray(&ptr2[idx+2*VecType::GetVecLength()]);
                scalarProduct+=res.horizontalSum();
                res=vec_ptr1.setFromArray(&ptr1[idx+3*VecType::GetVecLength()])*vec_ptr2.setFromArray(&ptr2[idx+3*VecType::GetVecLength()]);
                scalarProduct+=res.horizontalSum();
            }
        }

        for(unsigned long idx = nbValuesRounded4 ; idx < nbValues ; ++idx){
            scalarProduct += ptr1[idx] * ptr2[idx];
        }

        return scalarProduct;

    }


    inline RealType* ProductVecMat(RealType* vect, RealType* mat,
                                    const unsigned long nbRows, const unsigned long nbCols,
                                    const unsigned long nbValues){

        RealType* subMat = new RealType[nbValues];
        RealType* ptrRes = new RealType[nbRows];

        for(unsigned long idx = 0 ; idx < nbRows ; idx ++){
            memcpy(subMat, mat+(idx*nbCols), nbCols*sizeof(RealType));
            ptrRes[idx]=ScalarProduct(vect, subMat, nbValues);
        }

        return ptrRes;

    }

};




#endif // INABLAS_HPP
