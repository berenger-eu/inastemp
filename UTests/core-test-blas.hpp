///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef CORETESTBLAS_HPP
#define CORETESTBLAS_HPP

#include "InastempGlobal.h"
#include "UTester.hpp"
#include "Common/InaBlas.hpp"
#include "Common/InaUtils.hpp"

#include <cmath>
#include <cstring>
#include <cassert>

template < class VecType >
class TestBlas : public UTester< TestBlas< VecType > > {
    using Parent = UTester< TestBlas< VecType > >;

    using RealType = typename VecType::RealType;
    using MaskType = typename VecType::MaskType;

    void equalToVecType(const VecType vec,
                       const VecType inReal) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(vec.at(int(idx)), inReal.at(int(idx)));
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        RealType realsReals[VecType::GetVecLength()];
        inReal.storeInArray( realsReals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(reals[idx], realsReals[idx]);
        }

        alignas(512) RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        alignas(512) RealType realsRealsalign[VecType::GetVecLength()];
        inReal.storeInArray( realsRealsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(realsalign[idx], realsRealsalign[idx]);
        }
    }

    void equalToScalar(const VecType vec,
                       const RealType inReal) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(vec.at(int(idx)), inReal);
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(reals[idx], inReal);
        }

        alignas(512) RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(realsalign[idx], inReal);
        }
    }

    void equalToScalarMask(const MaskType vec,
                           const MaskType inReal) {
        VecType vecval = VecType::IfTrue(vec, VecType(1));
        VecType realval = VecType::IfTrue(inReal, VecType(1));
        equalToVecType(vecval, realval);
    }

    void equalToArray(const VecType vec,
                      const RealType inReals[]) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(vec.at(int(idx)), inReals[idx]);
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(reals[idx], inReals[idx]);
        }

        alignas(512) RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(realsalign[idx], inReals[idx]);
        }

        alignas(512) char reals_forcena_buffer[VecType::GetVecLength()*sizeof(RealType)+1];
        RealType* reals_forcena = reinterpret_cast<RealType*>(&reals_forcena_buffer[1]);
        vec.storeInArray( reals_forcena);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTEEQUAL(reals_forcena[idx], inReals[idx]);
        }
    }

    bool approxEqual(const float v1, const float v2) {
        return (std::abs(v1 - v2) / v2) <= 9.9999999999E-6f;
    }

    bool approxEqual(const double v1, const double v2) {
        return (std::abs(v1 - v2) / v2) <= 9.999999999999999E-12;
    }

    bool approxEqualLowAcc(const float v1, const float v2) {
        return (std::abs(v1 - v2) / v2) <= 9.9999999999E-2f;
    }

    bool approxEqualLowAcc(const double v1, const double v2) {
        return (std::abs(v1 - v2) / v2) <= 9.999999999999999E-5;
    }

    void approxEqualToScalar(const VecType vec,
                             const RealType inReal) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(vec.at(int(idx)), inReal));
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(reals[idx], inReal));
        }

        alignas(512) RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(realsalign[idx], inReal));
        }
    }

    void approxEqualToArray(const VecType vec,
                            const RealType inReals[]) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(vec.at(int(idx)), inReals[idx]));
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(reals[idx], inReals[idx]));
        }

        alignas(512) RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqual(realsalign[idx], inReals[idx]));
        }
    }

    void approxLowAccEqualToArray(const VecType vec,
                            const RealType inReals[]) {
        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqualLowAcc(vec.at(int(idx)), inReals[idx]));
        }

        RealType reals[VecType::GetVecLength()];
        vec.storeInArray( reals);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqualLowAcc(reals[idx], inReals[idx]));
        }

        alignas(512) RealType realsalign[VecType::GetVecLength()];
        vec.storeInAlignedArray( realsalign);

        for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
            UASSERTETRUE(approxEqualLowAcc(realsalign[idx], inReals[idx]));
        }
    }

    void TestBasic() {

        // ZERO
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(0);
            }

            alignas(512) char buff[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test = reinterpret_cast<RealType*>(&buff[1]);

            VecType vec(reals);

            InaBlas<VecType> blas{};

            blas.setZero(ptr_test, 37);

            equalToArray(vec, ptr_test);
            equalToScalar(VecType::GetZero(), ptr_test[36]);
        }


        // SCALAR
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(5);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            alignas(512) char buff[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test = reinterpret_cast<RealType*>(&buff[1]);

            blas.setScalar(ptr_test, 37, 5);

            equalToArray(vec, ptr_test);
            equalToScalar(VecType(5), ptr_test[36]);
        }

        // MULT BY SCALAR
        {
            alignas(512) RealType reals[VecType::GetVecLength()];
            for (size_t idx = 0; idx < size_t(VecType::GetVecLength()) ; ++idx) {
                reals[idx] = RealType(15);
            }

            VecType vec(reals);

            InaBlas<VecType> blas{};

            alignas(512) char buff[VecType::GetVecLength()*sizeof(RealType)+1];
            RealType* ptr_test = reinterpret_cast<RealType*>(&buff[1]);

            blas.setScalar(ptr_test, 37, 5);

            blas.VecMultScalar(ptr_test, 37, 3);

            equalToArray(vec, ptr_test);
            equalToScalar(VecType(15), ptr_test[36]);

        }

    }

    void SetTests() {
        Parent::AddTest(&TestBlas::TestBasic, "Basic test for blas functions");
    }

};

#endif
