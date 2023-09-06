#include "stdafx.h"
#include <wtypes.h>
#include <immintrin.h>
#include "avx_print.h"



void fir_avx512_fp16(int taps_samp, __bfloat16* barrayh, __m256i* startPtr, __m256h* filtered_data)
{
	__m512h x1_accumh = _mm512_setzero_ph();

	int count = 0;
	const int num_bits = 14;
	const int dc_offset = 32000;

	if (taps_samp > 0)
	{
		//taps_samp is always odd 
		int t;

		for (t = 0; t < taps_samp + 0; t += 2)
		{
			__m128i b1 = _mm_loadu_si16((barrayh + count++));

			//take lowest 16 bit value and broadcast - could do _mm512_mask_broadcastw_epi16 here for ZMM
			__m512h bh = _mm512_broadcastw_epi16(b1);

			__m128i b2 = _mm_loadu_si16((barrayh + count++));
			__m256h bh2 = _mm256_broadcastw_epi16(b2);

			bh = _mm512_inserti32x8(bh, bh2, 1);

			//load 32 channels of int16 RF data
			__m512i x1 = _mm512_loadu_epi16(((__m256i*)startPtr) + t);

			//shift off any zeros in LSBs
			x1 = _mm512_srli_epi16(x1, (short)(16 - num_bits));

			//remove any dc offset, critical for fp16 due to limited dynamic range
			x1 = _mm512_sub_epi16(x1, _mm512_set1_epi16(dc_offset >> (16 - num_bits)));
			//note:  value is now signed
			//printf("x1: ");	print256_i16(x1);

			//convert all 16 (signed if DC subtracted) channels to fp16
			__m512h x1hp = _mm512_cvtepi16_ph(x1);
			//__m256h x1hp = _mm256_cvtepu16_ph(x1);
			//printf("x1hp: "); print256_f16(x1hp);

			//multiply the fp16 RF data by the fp16 b filter coefficient and store into x1_accumh					
			x1_accumh = _mm512_fmadd_ph(x1hp, bh, x1_accumh);
		}

		//handle last (odd) tap and merge b1 and b2 accum 

		//merge b1 and b2 first
		__m256h accum1 = _mm512_extracti32x8_epi32(x1_accumh, 0);
		__m256h accum2 = _mm512_extracti32x8_epi32(x1_accumh, 1);
		__m256h accumCombined = _mm256_add_ph(accum1, accum2);

		__m128i b1 = _mm_loadu_si16((barrayh + count++));
		//take lowest 16 bit value and broadcast - could do _mm512_mask_broadcastw_epi16 here for ZMM
		__m256h bh = _mm256_broadcastw_epi16(b1);
		//load all 16 channels of int16 RF data
		__m256i x1 = _mm256_lddqu_si256(((__m256i*)startPtr) + t - 1);
		//shift off any zeros in LSBs
		x1 = _mm256_srli_epi16(x1, (short)(16 - num_bits));
		//remove any dc offset, critical for fp16 due to limited dynamic range
		x1 = _mm256_sub_epi16(x1, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));
		//note:  value is now signed
		//printf("x1: ");	print256_i16(x1);
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


void fir_avx256_fp16(int taps_samp, __bfloat16* barrayh, __m256i* startPtr, __m256h* filtered_data)
{

	__m256h x1_accumh = _mm256_setzero_ph();

	int count = 0;
	const int num_bits = 14;
	const int dc_offset = 32000;


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


//__attribute__((noinline))
void fir_x4_avx256_fp16(int taps_samp, __bfloat16* barrayh, __m256i* startPtr, __m256h* filtered_data, int num_bits, int dc_offset)
{
	__m256h x1_accumh = _mm256_setzero_ph();

	int count = 0;
	int total_taps = taps_samp + 1;

	for (int t = 0; t < total_taps; t++)
	{
		//when we computed taps_samp, we made sure that each group of 4 samples had the same number of taps, so this is safe
		//TODO:  optimize by repacking the barray to be continuous in memory
		

		unsigned short b1 = *((unsigned short*)(barrayh + count));
		unsigned short b2 = *((unsigned short*)(barrayh + count + 1 * total_taps));
		unsigned short b3 = *((unsigned short*)(barrayh + count + 2 * total_taps));
		unsigned short b4 = *((unsigned short*)(barrayh + count + 3 * total_taps));

		//print256_f16(*((__m256h*)(barrayh + count)));
		count++;

	
		__m256i bb2 = _mm256_set_epi16(b1, b1, b1, b1, b2, b2, b2, b2, b3, b3, b3, b3, b4, b4, b4, b4);

		//print256_f16(bb);

		//load 4 samples of 4 channels of int16 RF data

		//is this right?????  we are stepping 4 samples at a time within a single channel, which is wrong?
		__m256i x1 = _mm256_lddqu_si256((__m256i*)(((char*)startPtr) + t*8)); 

		x1 = _mm256_srli_epi16(x1, (short)(16 - num_bits));
		//remove any dc offset, critical for fp16 due to limited dynamic range
		x1 = _mm256_sub_epi16(x1, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));
		//note:  value is now signed

		//convert all 16 (signed if DC subtracted) channels to fp16
		__m256h x1hp = _mm256_cvtepi16_ph(x1);
		//printf("%d\n", t);
		//print256_f16(bb);
		//print256_f16(x1hp);

		//multiply the fp16 RF data by the fp16 b filter coefficient and store into x1_accumh					
		x1_accumh = _mm256_fmadd_ph(x1hp, _mm256_castsi256_ph(bb2), x1_accumh);
		
		//print256_f16(x1_accumh);
		//printf("\n");
	}
	//print256_f16(x1_accumh);
	_mm256_storeu_epi16((filtered_data), x1_accumh);

}

//optimized version that assumes the b coefficients are interleaved so a single _mm_loadu_ph can get 4
void fir_x4_avx256_fp16_opt(int taps_samp, __bfloat16* barrayh, __m256i* startPtr, __m256h* filtered_data, int num_bits, int dc_offset)
{
	__m256h x1_accumh = _mm256_setzero_ph();

	int count = 0;
	int total_taps = taps_samp + 1;

	for (int t = 0; t < total_taps; t++)
	{
		//load 4 packed coefficients, assumes that the precomputation has reorded them correctly
		__m128i bcoeffs = _mm_loadu_ph((barrayh + count));

		//shuffle the 4 values so that they fill a 16 element YMM regisgter
		__m256i idx = _mm256_set_epi16(3, 3, 3, 3, 2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0);
		__m256h b = _mm256_mask_permutex2var_epi16(_mm256_castsi128_si256(bcoeffs), 0xFFFFFFFF, idx, _mm256_castsi128_si256 (bcoeffs));

		count+=4;

		__m256i x1 = _mm256_lddqu_si256((__m256i*)(((char*)startPtr) + t * 8));

		x1 = _mm256_srli_epi16(x1, (short)(16 - num_bits));
		//remove any dc offset, critical for fp16 due to limited dynamic range
		x1 = _mm256_sub_epi16(x1, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));
		//note:  value is now signed

		//convert all 16 (signed if DC subtracted) channels to fp16
		__m256h x1hp = _mm256_cvtepi16_ph(x1);
	
		//multiply the fp16 RF data by the fp16 b filter coefficient and store into x1_accumh					
		x1_accumh = _mm256_fmadd_ph(x1hp, _mm256_castsi256_ph(b), x1_accumh);
	}

	_mm256_storeu_epi16((filtered_data), x1_accumh);

}