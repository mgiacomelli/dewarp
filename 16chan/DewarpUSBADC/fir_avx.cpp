#include "stdafx.h"
#include <wtypes.h>
#include <immintrin.h>


void fir_avx512_fp16(int taps_samp, __bfloat16* barrayh, __m256i* startPtr, __m256h* filtered_data, int num_bits, int dc_offset)
{
	__m512h x1_accumh = _mm512_setzero_ph();

	int count = 0;

	if (taps_samp > 0)
	{
		//taps_samp is always odd 
		int t;

		for (t = 0; t < taps_samp + 0; t += 2)
		{
			__m128i b1 = _mm_loadu_si16((barrayh + count++));

			//take lowest 16 bit value and broadcast 
			__m512h bh = _mm512_broadcastw_epi16(b1);

			__m128i b2 = _mm_loadu_si16((barrayh + count++));
			__m256h bh2 = _mm256_broadcastw_epi16(b2);

			bh = _mm512_inserti32x8(bh, bh2, 1);

			//load 32 channels of int16 RF data
			__m512i x1 = _mm512_loadu_epi16(startPtr + t);

			//shift off any zeros in LSBs
			x1 = _mm512_srli_epi16(x1, (short)(16 - num_bits));

			//remove any dc offset, critical for fp16 due to limited dynamic range
			x1 = _mm512_sub_epi16(x1, _mm512_set1_epi16(dc_offset >> (16 - num_bits)));
			//note:  value is now signed

			//convert all 16 (signed if DC subtracted) channels to fp16
			__m512h x1hp = _mm512_cvtepi16_ph(x1);
			//__m256h x1hp = _mm256_cvtepu16_ph(x1);

			//multiply the fp16 RF data by the fp16 b filter coefficient and store into x1_accumh					
			x1_accumh = _mm512_fmadd_ph(x1hp, bh, x1_accumh);
		}

		//handle last (odd) tap and merge b1 and b2 accum 

		//merge b1 and b2 first
		__m256h accum1 = _mm512_extracti32x8_epi32(x1_accumh, 0);
		__m256h accum2 = _mm512_extracti32x8_epi32(x1_accumh, 1);
		__m256h accumCombined = _mm256_add_ph(accum1, accum2);

		__m128i b1 = _mm_loadu_si16((barrayh + count++));
		//take lowest 16 bit value and broadcast 
		__m256h bh = _mm256_broadcastw_epi16(b1);
		//load all 16 channels of int16 RF data
		__m256i x1 = _mm256_lddqu_si256(((__m256i*)startPtr) + t - 1);
		//shift off any zeros in LSBs
		x1 = _mm256_srli_epi16(x1, (short)(16 - num_bits));
		//remove any dc offset, critical for fp16 due to limited dynamic range
		x1 = _mm256_sub_epi16(x1, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));
		//note:  value is now signed

		//convert all 16 (signed if DC subtracted) channels to fp16
		__m256h x1hp = _mm256_cvtepi16_ph(x1);

		accumCombined = _mm256_fmadd_ph(x1hp, bh, accumCombined);
		_mm256_storeu_epi16((short*)(filtered_data), accumCombined);
	}
	else
	{
		//optimize case where we don't really filter (1 tap filter)
		count++;
		__m256i x1 = _mm256_lddqu_si256(((__m256i*)startPtr));
		x1 = _mm256_srli_epi16(x1, (short)(16 - num_bits));
		//remove any dc offset, critical for fp16 due to limited dynamic range
		x1 = _mm256_sub_epi16(x1, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));

		//convert to fp16 and store directly

		_mm256_storeu_epi16((short*)(filtered_data), _mm256_cvtepi16_ph(x1));
	}
}

void fir_avx256_fp16(int taps_samp, __bfloat16* barrayh, __m256i* startPtr, __m256h* filtered_data, int num_bits, int dc_offset)
{

	__m256h x1_accumh = _mm256_setzero_ph();

	int count = 0;

	int t;
	for (t = 0; t < taps_samp + 1; t++)
	{
		__m128i b1 = _mm_loadu_si16((barrayh + count++));
		//take lowest 16 bit value and broadcast 
		__m256h bh = _mm256_broadcastw_epi16(b1);
		//load all 16 channels of int16 RF data
		__m256i x1 = _mm256_lddqu_si256(((__m256i*)startPtr) + t);
		//shift off any zeros in LSBs
		x1 = _mm256_srli_epi16(x1, (short)(16 - num_bits));
		//remove any dc offset, critical for fp16 due to limited dynamic range
		x1 = _mm256_sub_epi16(x1, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));
		//note:  value is now signed

		//convert all 16 (signed if DC subtracted) channels to fp16
		__m256h x1hp = _mm256_cvtepi16_ph(x1);
		//__m256h x1hp = _mm256_cvtepu16_ph(x1);

		//multiply the fp16 RF data by the fp16 b filter coefficient and store into x1_accumh					
		x1_accumh = _mm256_fmadd_ph(x1hp, bh, x1_accumh);
	}

	_mm256_storeu_epi16((filtered_data), x1_accumh);
}