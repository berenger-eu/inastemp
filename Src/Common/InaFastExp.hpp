///////////////////////////////////////////////////////////////////////////
// Inastemp - Berenger Bramas MPCDF - 2016
// Under MIT Licence, please you must read the LICENCE file.
///////////////////////////////////////////////////////////////////////////
#ifndef INAFASTEXP_HPP
#define INAFASTEXP_HPP

#include <array>
#include <cmath>

/**
 * This class provides the methods and constant to
 * compute the fast exp functions.
 * From the paper: Fast Exponential Computation on SIMD Architectures.
 * The constant from Remez polynomials are generated
 * with remez.sce (scilab) see documentation to know more.
 */
class InaFastExp {
    static const long int S64 = (1L << 52);
    static const long int S32 = (1L << 23);

public:
    constexpr static double CoeffLog2E() {
        // const double Euler = 2.71828182845904523536028747135266249775724709369995;
        // return std::log2(Euler());
        return 1.442695040888963407359924681001892137426645954153;
    }

    constexpr static double CoeffA64() {
        return double(S64);
    }

    constexpr static double CoeffB64() {
        return double(S64 * 1023);
    }

    constexpr static double CoeffA32() {
        return double(S32);
    }

    constexpr static double CoeffB32() {
        return double(S32 * 127);
    }

    constexpr static std::array< double, 3 > GetCoefficient3() {
        return { { -2.47142816509723700288e-03,
                   3.48943001636461247461e-01,
                   -3.44000145306266491563e-01 } };
    }

    constexpr static std::array< double, 4 > GetCoefficient4() {
        return { { 1.06906116358144185133e-04,
                   3.03543677780836240743e-01,
                   -2.24339532327269441936e-01,
                   -7.92041454535668681958e-02 } };
    }

    constexpr static std::array< double, 5 > GetCoefficient5() {
        return { { -3.70138142771437266806e-06,
                   3.07033820309224325662e-01,
                   -2.41638288055762540107e-01,
                   -5.16904731562965388814e-02,
                   -1.36976563343097993558e-02 } };
    }

    constexpr static std::array< double, 6 > GetCoefficient6() {
        return { { 1.06823753710239477000e-07,
                   3.06845249656632845792e-01,
                   -2.40139721982230797126e-01,
                   -5.58662282412822480682e-02,
                   -8.94283890931273951763e-03,
                   -1.89646052380707734290e-03 } };
    }

    constexpr static std::array< double, 7 > GetCoefficient7() {
        return { { -2.64303273610414963822e-09,
                   3.06853075372807815313e-01,
                   -2.40230549677691723742e-01,
                   -5.54802224547989303316e-02,
                   -9.68497459444197204836e-03,
                   -1.23843111224273085859e-03,
                   -2.18892247566917477666e-04 } };
    }

    constexpr static std::array< double, 8 > GetCoefficient8() {
        return { { 5.72265234348656066133e-11,
                   3.06852812183173784266e-01,
                   -2.40226356058427820139e-01,
                   -5.55053022725605083032e-02,
                   -9.61350625581030605871e-03,
                   -1.34302437845634000529e-03,
                   -1.42962470418959216190e-04,
                   -2.16607474999407558923e-05 } };
    }

    constexpr static std::array< double, 9 > GetCoefficient9() {
        return { { -1.10150186041739869460e-12,
                   3.06852819617161765020e-01,
                   -2.40226511645233870018e-01,
                   -5.55040609720754696266e-02,
                   -9.61837182960864275905e-03,
                   -1.33266405715993271723e-03,
                   -1.55186852765468104613e-04,
                   -1.41484352491262699514e-05,
                   -1.87582286605066256753e-06 } };
    }

    constexpr static std::array< double, 10 > GetCoefficient10() {
        return { { 1.79040320992239805871e-11,
                   3.06852815055756844576e-01,
                   -2.40226385506041861806e-01,
                   -5.55053584940081654042e-02,
                   -9.61174262279892825667e-03,
                   -1.35164210003994454852e-03,
                   -1.23291147980286128769e-04,
                   -4.53940620364305641833e-05,
                   1.46363500589519947862e-05,
                   -3.63750326480946818984e-06 } };
    }

    constexpr static std::array< double, 11 > GetCoefficient11() {
        return { { 7.32388148129676088418e-13,
                   3.06852819216552274995e-01,
                   -2.40226499275945526435e-01,
                   -5.55042073859920090384e-02,
                   -9.61749102796678571881e-03,
                   -1.33571753728242812072e-03,
                   -1.48718480159015542822e-04,
                   -2.26598047213222231406e-05,
                   4.91492761180322572151e-06,
                   -3.00875847392884227107e-06,
                   5.68126156224525271282e-07 } };
    }
};

#endif
