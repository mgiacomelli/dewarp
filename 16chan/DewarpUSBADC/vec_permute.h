#pragma once


//#define DEBUG
//16 way interleaved data, 16 samples at a time into 16 individual buffers, 32 total permutes (2 per output sample)
void deinterleave16x16(__m256i* input, UINT16* output, int outstride)
{
	/*deinterleave 16 YMM registers containing 16, 16 bit samples
	 Proceed in 4 passes, first creates pairs of sequential samples,
	 then second combines 2 pairs to make quads, then third combines
	 2 quads to make an octal.  Finally 4th passe puts 16 samples together
	*/

	__m512i intermediatearray[8];

	__m512i intermediatearray2a[4];
	__m512i intermediatearray2b[4];

	__m512i intermediatearray3a[2];
	__m512i intermediatearray3b[2];
	__m512i intermediatearray3c[2];
	__m512i intermediatearray3d[2];


	//__m256h tmp = _mm256_set_ph(240, 224, 208, 192, 176, 160, 144, 128, 112, 96, 80, 64, 48, 32, 16, 0);
#ifdef DEBUG
	printf("\n--------16 Way Interleaved DATA-----------------\n");

	for (int i = 0; i < 16; i++)
	{
		//    tosavearray[i] = _mm256_add_ph(tmp, _mm256_set1_ph((float)i));
		print256_i16(input[i]);
	}


	printf("\n--------First Round Merge-----------------\n");
#endif
	//go pairwise and combine sequential samples
	__m512i idx1 = _mm512_set_epi16(47, 15, 46, 14, 45, 13, 44, 12, 43, 11, 42, 10, 41, 9, 40, 8, 39, 7, 38, 6, 37, 5, 36, 4, 35, 3, 34, 2, 33, 1, 32, 0);
	for (int i = 0; i < 16; i += 2)
	{
		//todo:  see if can replace with _mm512_inserti32x8 which has double the throughput of permute
		__m512i tmp = _mm512_mask_permutex2var_epi16(_mm512_castsi256_si512(input[i]), 0xFFFFFFFF, idx1, _mm512_castsi256_si512(input[i + 1]));
		intermediatearray[i / 2] = tmp;
#ifdef DEBUG
		print512_i16(tmp);
#endif

	}

#ifdef DEBUG
	printf("\n--------Second Round Merge-----------------\n");
#endif
	__m512i idx2 = _mm512_set_epi16(47, 46, 15, 14, 45, 44, 13, 12, 43, 42, 11, 10, 41, 40, 9, 8, 39, 38, 7, 6, 37, 36, 5, 4, 35, 34, 3, 2, 33, 32, 1, 0);
	// +16 to get the second half
	__m512i idx3 = _mm512_set_epi16(63, 62, 31, 30, 61, 60, 29, 28, 59, 58, 27, 26, 57, 56, 25, 24, 55, 54, 23, 22, 53, 52, 21, 20, 51, 50, 19, 18, 49, 48, 17, 16);

	for (int i = 0; i < 8; i += 2)
	{
		__m512i tmp = _mm512_mask_permutex2var_epi16(intermediatearray[i], 0xFFFFFFFF, idx2, intermediatearray[i + 1]);
		intermediatearray2a[i / 2] = tmp;
#ifdef DEBUG
		print512_i16(tmp);
#endif

		__m512i tmp2 = _mm512_mask_permutex2var_epi16(intermediatearray[i], 0xFFFFFFFF, idx3, intermediatearray[i + 1]);
		intermediatearray2b[i / 2] = tmp2;
#ifdef DEBUG
		print512_i16(tmp2);
#endif

	}

#ifdef DEBUG
	printf("\n--------Third Round Merge-----------------\n");
#endif

	//__m512i idx4 = _mm512_set_epi16(47, 46, 45, 15, 14, 13, 12, 44, 43, 42, 41, 40, 11, 10, 9, 8, 39, 38, 37, 36, 7, 6, 5, 4, 35, 34, 33, 32, 3, 2, 1, 0);
	//__m512i idx5 = _mm512_set_epi16(63, 62, 61, 31, 30, 29, 28, 60, 59, 58, 57, 56, 27, 26, 25, 24, 55, 54, 53, 52, 23, 22, 21, 20, 51, 50, 49, 48, 19, 18, 17, 16);
	__m512i idx4 = _mm512_set_epi16(47, 46, 45, 44, 15, 14, 13, 12, 43, 42, 41, 40, 11, 10, 9, 8, 39, 38, 37, 36, 7, 6, 5, 4, 35, 34, 33, 32, 3, 2, 1, 0);
	__m512i idx5 = _mm512_set_epi16(63, 62, 61, 60, 31, 30, 29, 28, 59, 58, 57, 56, 27, 26, 25, 24, 55, 54, 53, 52, 23, 22, 21, 20, 51, 50, 49, 48, 19, 18, 17, 16);

	for (int i = 0; i < 4; i += 2)
	{
		__m512i tmp = _mm512_mask_permutex2var_epi16(intermediatearray2a[i], 0xFFFFFFFF, idx4, intermediatearray2a[i + 1]);
		intermediatearray3a[i / 2] = tmp;
#ifdef DEBUG
		print512_i16(tmp);
#endif

		__m512i tmp2 = _mm512_mask_permutex2var_epi16(intermediatearray2a[i], 0xFFFFFFFF, idx5, intermediatearray2a[i + 1]);
		intermediatearray3b[i / 2] = tmp2;
#ifdef DEBUG
		print512_i16(tmp2);
#endif
	}

	for (int i = 0; i < 4; i += 2)
	{
		__m512i tmp = _mm512_mask_permutex2var_epi16(intermediatearray2b[i], 0xFFFFFFFF, idx4, intermediatearray2b[i + 1]);

		intermediatearray3c[i / 2] = tmp;

		__m512i tmp2 = _mm512_mask_permutex2var_epi16(intermediatearray2b[i], 0xFFFFFFFF, idx5, intermediatearray2b[i + 1]);

		intermediatearray3d[i / 2] = tmp2;
#ifdef DEBUG
		print512_i16(tmp);
		print512_i16(tmp2);
#endif
	}

#ifdef DEBUG
	//printf("\n--------Last Round Merge-----------------\n [===============================  Sample 0  ==================================================] [===============================  Sample 1  ==================================================]\n");
	printf("\n--------Last Round Merge-----------------\n [===============================  Channel N  =================================================]\n");
#endif
	__m512i idx6 = _mm512_set_epi16(47, 46, 45, 44, 43, 42, 41, 40, 15, 14, 13, 12, 11, 10, 9, 8, 39, 38, 37, 36, 35, 34, 33, 32, 7, 6, 5, 4, 3, 2, 1, 0);
	__m512i idx7 = _mm512_set_epi16(63, 62, 61, 60, 59, 58, 57, 56, 31, 30, 29, 28, 27, 26, 25, 24, 55, 54, 53, 52, 51, 50, 49, 48, 23, 22, 21, 20, 19, 18, 17, 16);

	{
		__m512i tmp = _mm512_mask_permutex2var_epi16(intermediatearray3a[0], 0xFFFFFFFF, idx6, intermediatearray3a[1]);
		//output[0] = tmp;
		*(__m256i*)(output) =  _mm512_extracti32x8_epi32(tmp, 0);
		*(__m256i*)(output+outstride) =  _mm512_extracti32x8_epi32(tmp, 1);

		tmp = _mm512_mask_permutex2var_epi16(intermediatearray3a[0], 0xFFFFFFFF, idx7, intermediatearray3a[1]);
		//output[1] = tmp;
		*(__m256i*)(output + outstride*2) =  _mm512_extracti32x8_epi32(tmp, 0);

		//THIS VALUE IS WRONG (CONTENTS SCRAMBLED)
		*(__m256i*)(output + outstride*3) =  _mm512_extracti32x8_epi32(tmp, 1);


		__m512i tmp2 = _mm512_mask_permutex2var_epi16(intermediatearray3b[0], 0xFFFFFFFF, idx6, intermediatearray3b[1]);
		//output[2] = tmp2;
		*(__m256i*)(output + outstride * 4) = _mm512_extracti32x8_epi32(tmp2, 0);
		*(__m256i*)(output + outstride * 5) = _mm512_extracti32x8_epi32(tmp2, 1);

		tmp2 = _mm512_mask_permutex2var_epi16(intermediatearray3b[0], 0xFFFFFFFF, idx7, intermediatearray3b[1]);
		//output[3] = tmp2;
		*(__m256i*)(output + outstride * 6) = _mm512_extracti32x8_epi32(tmp2, 0);
		*(__m256i*)(output + outstride * 7) = _mm512_extracti32x8_epi32(tmp2, 1);

		__m512i tmp3 = _mm512_mask_permutex2var_epi16(intermediatearray3c[0], 0xFFFFFFFF, idx6, intermediatearray3c[1]);
		//output[4] = tmp3;
		*(__m256i*)(output + outstride * 8) = _mm512_extracti32x8_epi32(tmp3, 0);
		*(__m256i*)(output + outstride * 9) = _mm512_extracti32x8_epi32(tmp3, 1);

		tmp3 = _mm512_mask_permutex2var_epi16(intermediatearray3c[0], 0xFFFFFFFF, idx7, intermediatearray3c[1]);
		//output[5] = tmp3;
		*(__m256i*)(output + outstride * 10) = _mm512_extracti32x8_epi32(tmp3, 0);
		*(__m256i*)(output + outstride * 11) = _mm512_extracti32x8_epi32(tmp3, 1);


		__m512i tmp4 = _mm512_mask_permutex2var_epi16(intermediatearray3d[0], 0xFFFFFFFF, idx6, intermediatearray3d[1]);
		//output[6] = tmp4;
		*(__m256i*)(output + outstride * 12) = _mm512_extracti32x8_epi32(tmp4, 0);
		*(__m256i*)(output + outstride * 13) = _mm512_extracti32x8_epi32(tmp4, 1);

		tmp4 = _mm512_mask_permutex2var_epi16(intermediatearray3d[0], 0xFFFFFFFF, idx7, intermediatearray3d[1]);
		//output[7] = tmp4;
		*(__m256i*)(output + outstride * 14) = _mm512_extracti32x8_epi32(tmp4, 0);
		*(__m256i*)(output + outstride * 15) = _mm512_extracti32x8_epi32(tmp4, 1);

	}

#ifdef DEBUG
	for (int i = 0; i < 16; i++)
	{
		print256_i16(*(__m256i*)(output + i*outstride));
	}
#endif
}
