#pragma once
#include "stdafx.h"
#include <wtypes.h>
#include <immintrin.h>



union U512 {
    __m512i vi;
    __m512 v;
    INT16 i16[32];
    INT32 i32[16];
    float f32[16];
};

union U256f {
    __m256 v;
    float a[8];
};

union U256i {
    __m256i v;
    INT16 a[16];
    INT32 b[8];
};

void print256_f32(const __m256 v)
{
    const U256f u = { v };

    for (int i = 0; i < 8; ++i)
        printf("%f\n", u.a[i]);
}

void print256_f16(const __m256h v)
{
    //printf doesn't support fp16, so convert ot fp32
    __m512 mm512_fp32 = _mm512_cvtph_ps((__m256i)v);

    U512 u = { };
    u.v = mm512_fp32;

    for (int i = 0; i < 16; ++i)
    {
        printf("% 05f ", u.f32[i]);
    }
    printf("\n");

}

void print512_f32(const __m512 v)
{
    U512 u = {  };
    u.v = v;

    for (int i = 0; i < 16; ++i)
        printf("% 05f  ", u.f32[i]);
}

void print512_f16(const __m512i v)
{
    __m256i a = _mm512_extracti64x4_epi64(v, 0);
    __m512 mm512_fp32 = _mm512_cvtph_ps(a);
    U512 u1 = {  };
    u1.v = mm512_fp32;

    __m256i b = _mm512_extracti64x4_epi64(v, 1);
    mm512_fp32 = _mm512_cvtph_ps(b);
    U512 u2 = {  };
    u2.v = mm512_fp32;


    for (int i = 0; i < 16; ++i)
        printf("% 05f  ", u1.f32[i]);
    for (int i = 0; i < 16; ++i)
        printf("% 05f  ", u2.f32[i]);
    printf("\n");
}


void print128_f16(const __m128h v)
{
    //printf doesn't support fp16, so convert ot fp32
    __m256 fp32 = _mm256_cvtph_ps(v);

    const U256f u = { fp32 };

    for (int i = 0; i < 8; ++i)
    {
        printf("% 05f  ", u.a[i]);
    }
    printf("\n");
}

void print512_i16(__m512i v)
{

    U512 u = { };
    u.vi = v;

    for (int i = 0; i < 32; ++i)
        printf("% 05hd ", u.i16[i]);
    printf("\n");

}

void print256_i16(__m256i v)
{
    const U256i u = { v };

    for (int i = 0; i < 16; ++i)
        printf("% 05hd ", u.a[i]);
    printf("\n");
}
void print256_u16(__m256i v)
{
    const U256i u = { v };

    for (int i = 0; i < 16; ++i)
        printf("% 05hu ", u.a[i]);
    printf("\n");
}

void print256i_i32(__m256i v)
{
    const U256i u = { v };

    for (int i = 0; i < 8; ++i)
        printf("% 05d ", u.b[i]);
    printf("\n");
}

void print128_i16(__m128i packed1)
{
    printf("%i\n", _mm_extract_epi16(packed1, 0));
    printf("%i\n", _mm_extract_epi16(packed1, 1));
    printf("%i\n", _mm_extract_epi16(packed1, 2));
    printf("%i\n", _mm_extract_epi16(packed1, 3));
    printf("%i\n", _mm_extract_epi16(packed1, 4));
    printf("%i\n", _mm_extract_epi16(packed1, 5));
    printf("%i\n", _mm_extract_epi16(packed1, 6));
    printf("%i\n", _mm_extract_epi16(packed1, 7));
}
