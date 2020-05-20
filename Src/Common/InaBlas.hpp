///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAVECINTERFACE_HPP
#define INAVECINTERFACE_HPP

#include <cstdio>
#include <vector>
#include <cstring>

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


    inline RealType* ProductVecMat(RealType* vect, RealType** mat,
        const unsigned long nbRows,
        const unsigned long nbValues){

        RealType* ptrRes = new RealType[nbRows];

        for(unsigned long idx = 0 ; idx < nbRows ; idx ++){
            ptrRes[idx]=ScalarProduct(vect, mat[idx], nbValues);
        }

        return ptrRes;

    }


    inline RealType** ProductMatMat(RealType** mat1, RealType** mat2,
        const unsigned long nbRows1, const unsigned long nbCols,
        const unsigned long nbCols2){

        //RealType* subMat1 = new RealType[nbCols];
        RealType* subMat2 = new RealType[nbCols];

        RealType **matRes = new RealType*[nbRows1];
        for (unsigned long i=0; i < nbRows1; i++)
        matRes[i] = new RealType[nbCols2];

        for(unsigned long idx_row = 0 ; idx_row < nbRows1 ; idx_row ++){
            for(unsigned long idx_col = 0 ; idx_col < nbCols2 ; idx_col ++){

                // TODO: vec submat2
                for(unsigned long idx_row2 = 0 ; idx_row2 < nbCols ; idx_row2 ++){
                    subMat2[idx_row2] = mat2[idx_row2][idx_col];
                }

                matRes[idx_row][idx_col]=ScalarProduct(mat1[idx_row], subMat2, nbCols);
            }
        }

        return matRes;

    }


    inline RealType* ProductMatMat2(RealType* mat1, RealType* mat2,
                                    const unsigned long nbRows1, const unsigned long nbCols,
                                    const unsigned long nbCols2){

        RealType* matRes = new RealType[nbRows1*nbCols2];
        RealType* subMat1 = new RealType[nbCols];
        RealType* subMat2 = new RealType[nbCols];

        // transposee to simplify the multiplication
        RealType* matTrans = TransposeeOpti(mat2, nbCols, nbCols2);

        for(unsigned long idx=0 ; idx < nbRows1 ; idx++){
            memcpy(subMat1, mat1+(idx*nbCols), nbCols*sizeof(RealType));
            for(unsigned long idx2=0 ; idx2 < nbCols2 ; idx2++){
                memcpy(subMat2, matTrans+(idx2*nbCols), nbCols*sizeof(RealType));
                matRes[idx*nbCols2+idx2]=ScalarProduct(subMat1, subMat2, nbCols);
            }
        }

        return matRes;

    }


    inline RealType* Transposee(RealType* mat, const unsigned long nbRows, const unsigned long nbCols){

        VecType vecMat;
        RealType* matRes = new RealType[nbRows*nbCols];

        const unsigned long nbRowsTr = nbCols;
        const unsigned long nbColsTr = nbRows;

        int indirect[VecType::GetVecLength()];

        const unsigned long nbValuesRounded = nbRows - (nbRows%VecType::GetVecLength());

        for(unsigned long idxCol = 0 ; idxCol < nbCols ; idxCol++) {
            for(unsigned long idxRow = 0 ; idxRow < nbValuesRounded ; idxRow += VecType::GetVecLength()){
                for(unsigned long idx = 0; idx<size_t(VecType::GetVecLength()); idx++){
                    indirect[idx] = static_cast<int>((idx+idxRow)*nbCols+idxCol);
                }
                vecMat.setFromIndirectArray(mat, indirect);
                vecMat.storeInArray(&matRes[idxCol*nbColsTr+idxRow]);
            }
        }

        for(unsigned long idxCol = 0 ; idxCol < nbRowsTr ; idxCol++) {
            for(unsigned long idxRow = nbValuesRounded; idxRow < nbRows ; idxRow++){
                matRes[idxCol*nbColsTr+idxRow] = mat[idxRow*nbRowsTr+idxCol];
            }
        }

        return matRes;

    }


    inline RealType* TransposeeOpti(RealType* mat, const unsigned long nbRows, const unsigned long nbCols){

        VecType vecMat, vecMatTr;
        RealType* matRes = new RealType[nbRows*nbCols];

        const unsigned long nbRowsTr = nbCols;
        const unsigned long nbColsTr = nbRows;

        int indirect[4*VecType::GetVecLength()];

        const unsigned long nbValuesRounded4 = nbRows - (nbRows%(4*VecType::GetVecLength()));

        for(unsigned long idxCol = 0 ; idxCol < nbCols ; idxCol++) {
            for(unsigned long idxRow = 0 ; idxRow < nbValuesRounded4 ; idxRow += 4*VecType::GetVecLength()){
                for(unsigned long idx = 0; idx<size_t(4*VecType::GetVecLength()); idx++){
                    indirect[idx] = static_cast<int>((idx+idxRow)*nbCols+idxCol);
                }
                vecMat.setFromIndirectArray(mat, &indirect[0]);
                vecMat.storeInArray(&matRes[idxCol*nbColsTr+idxRow]);
                vecMat.setFromIndirectArray(mat, &indirect[VecType::GetVecLength()]);
                vecMat.storeInArray(&matRes[idxCol*nbColsTr+idxRow+VecType::GetVecLength()]);
                vecMat.setFromIndirectArray(mat, &indirect[2*VecType::GetVecLength()]);
                vecMat.storeInArray(&matRes[idxCol*nbColsTr+idxRow+2*VecType::GetVecLength()]);
                vecMat.setFromIndirectArray(mat, &indirect[3*VecType::GetVecLength()]);
                vecMat.storeInArray(&matRes[idxCol*nbColsTr+idxRow+3*VecType::GetVecLength()]);
            }
        }

        for(unsigned long idxCol = 0 ; idxCol < nbRowsTr ; idxCol++) {
            for(unsigned long idxRow = nbValuesRounded4; idxRow < nbRows ; idxRow++){
                matRes[idxCol*nbColsTr+idxRow] = mat[idxRow*nbRowsTr+idxCol];
            }
        }

        return matRes;

    }

};



#endif // INABLAS_HPP
