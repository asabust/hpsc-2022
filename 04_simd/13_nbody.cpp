#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <immintrin.h>

int main()
{
  const int N = 8;
  float x[N], y[N], m[N], fx[N], fy[N], idx[N];
  for (int i = 0; i < N; i++)
  {
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
    idx[i] = i;
  }
  __m256 xvec = _mm256_load_ps(x);
  __m256 yvec = _mm256_load_ps(y);
  __m256 mvec = _mm256_load_ps(m);
  __m256 idxvec = _mm256_load_ps(idx);

  for (int i = 0; i < N; i++)
  {

    __m256 fxvec = _mm256_set1_ps(0);
    __m256 fyvec = _mm256_set1_ps(0);
    __m256 xivec = _mm256_set1_ps(x[i]);
    __m256 yivec = _mm256_set1_ps(y[i]);

    __m256 rx = _mm256_sub_ps(xivec, xvec);
    __m256 ry = _mm256_sub_ps(yivec, yvec);
    __m256 r = _mm256_rsqrt_ps(
        _mm256_add_ps(
            _mm256_mul_ps(rx, rx),
            _mm256_mul_ps(ry, ry)));
    __m256 temp = _mm256_mul_ps(_mm256_mul_ps(r, r), _mm256_mul_ps(r, mvec));
    fxvec = _mm256_sub_ps(fxvec, _mm256_mul_ps(rx, temp));
    fyvec = _mm256_sub_ps(fyvec, _mm256_mul_ps(ry, temp));

    __m256 ivec = _mm256_set1_ps(i);
    __m256 mask = _mm256_cmp_ps(idxvec, ivec, _CMP_EQ_OQ);
    fxvec = _mm256_blendv_ps(fxvec, _mm256_setzero_ps(), mask);
    fyvec = _mm256_blendv_ps(fxvec, _mm256_setzero_ps(), mask);
    _mm256_store_ps(fx, fxvec);
    _mm256_store_ps(fy, fyvec);
    // for (int j = 0; j < N; j++)
    // {
    //   if (i != j)
    //   {
    //     float rx = x[i] - x[j];
    //     float ry = y[i] - y[j];
    //     float r = std::sqrt(rx * rx + ry * ry);
    //     fx[i] = f[i] -  rx * m[j] / (r * r * r);
    //     fy[i] -= ry * m[j] / (r * r * r);
    //   }
    //}

    printf("%d %g %g\n", i, fx[i], fy[i]);
  }
}
