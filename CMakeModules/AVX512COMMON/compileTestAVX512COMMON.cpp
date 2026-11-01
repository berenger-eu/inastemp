#include <immintrin.h>

int main() {
    const __m512d valuesDouble = _mm512_set1_pd(-1.5);
    const __m512 valuesFloat = _mm512_set1_ps(-1.5f);
    const __m512d roundedDouble = _mm512_roundscale_pd(valuesDouble, _MM_FROUND_TO_NEG_INF);
    const __m512 roundedFloat = _mm512_roundscale_ps(valuesFloat, _MM_FROUND_TO_NEG_INF);
    return _mm512_cmp_pd_mask(roundedDouble, _mm512_setzero_pd(), _CMP_EQ_OQ)
         + _mm512_cmp_ps_mask(roundedFloat, _mm512_setzero_ps(), _CMP_EQ_OQ);
}
