#pragma once

//#define DEBUG 

//permute 4 pixels (0->3) of each 
// [A4B4C4D4 A3B3C3D3 A2B2C2D2 A1B1C1D1]_0 
// [A4B4C4D4 A2B2C2D2 A1B1C1D1 A3B3C3D3 ]_1
// ...
// into
//[A1B1C1D1_0 A1B1C1D1_1 A1B1C1D1_2 A1B1C1D1_3]
//[A2B2C2D2_0 A2B2C2D2_1 A2B2C2D2_2 A2B2C2D2_3]
//[A3B3C3D3_0 A3B3C3D3_1 A3B3C3D3_2 A3B3C3D3_3]
//[A4B4C4D4_0 A4B4C4D4_1 A4B4C4D4_2 A4B4C4D4_3]

//note that the input samples are in reverse order (3->0) but are returned in order (0->3)

void deinterleave4x4(__m256i* input, __m256i* output)
{
#ifdef DEBUG
	printf("\n--------4 Way Interleaved DATA, 4 pixels of 4 taps each-----------------\n");
	printf("\n");
	printf("         [ ===========================================  YMM  =========================================   ] \n");
	printf("         [   TAP 4                   TAP 3                   TAP 2                   TAP 1	               ] \n");
	printf("         [   A     B     C     D     A     B     C     D     A     B     C     D     A     B     C     D ] \n");

	for (int i = 0; i < 4; i++)
	{
		printf("Sample %d: ", i);
		print256_i16(input[i]);
	}

	printf("\n--------First Round Merge into ZMM-----------------\n");
#endif

	__m512i idx = _mm512_set_epi16(47, 46, 45, 44, 15, 14, 13, 12, 43, 42, 41, 40, 11, 10, 9, 8, 39, 38, 37, 36, 7, 6, 5, 4, 35, 34, 33, 32, 3, 2, 1, 0);
	__m512i AB512 = _mm512_mask_permutex2var_epi16(_mm512_castsi256_si512(input[0]), 0xFFFFFFFF, idx, _mm512_castsi256_si512(input[1]));
	__m512i CD512 = _mm512_mask_permutex2var_epi16(_mm512_castsi256_si512(input[2]), 0xFFFFFFFF, idx, _mm512_castsi256_si512(input[3]));
	//__m512i AB512 = _mm512_inserti32x8(_mm512_castsi256_si512(input[0]), input[1], 1);
	//__m512i CD512 = _mm512_inserti32x8(_mm512_castsi256_si512(input[2]), input[3], 1);

#ifdef DEBUG
	print512_i16(AB512);
	print512_i16(CD512);

	printf("\n--------Resulting DATA in ZMM-----------------\n");
#endif

	 idx = _mm512_set_epi16(47, 46, 45, 44, 43, 42, 41, 40, 15, 14, 13, 12, 11, 10, 9, 8, 39, 38, 37, 36, 35, 34, 33, 32, 7, 6, 5, 4, 3, 2, 1, 0);

	__m512i tmp = _mm512_mask_permutex2var_epi16(AB512, 0xFFFFFFFF, idx, CD512);

#ifdef DEBUG
	print512_i16(tmp);
#endif

	output[3] = _mm512_extracti32x8_epi32(tmp, 0);
	output[2] = _mm512_extracti32x8_epi32(tmp, 1);

	idx = _mm512_set_epi16(63, 62, 61, 60, 59, 58, 57, 56, 31, 30, 29, 28, 27, 26, 25, 24, 55, 54, 53, 52, 51, 50, 49, 48, 23, 22, 21, 20, 19, 18, 17, 16);
	tmp = _mm512_mask_permutex2var_epi16(AB512, 0xFFFFFFFF, idx, CD512);

#ifdef DEBUG
	print512_i16(tmp);
#endif

	output[1] = _mm512_extracti32x8_epi32(tmp, 0);
	output[0] = _mm512_extracti32x8_epi32(tmp, 1);




#ifdef DEBUG

	printf("\n");
	printf("         [ =============================================  YMM  ========================================= ] \n");
	printf("         [   Sample 1                Sample 2                Sample 3                Sample 4            ] \n");
	printf("         [   A     B     C     D     A     B     C     D     A     B     C     D     A     B     C     D ] \n");
	//printf("      [   A1    B1    C1    D1    A2    B2    C2    D2    A3    B3    C3    D3    A4    B4    C4    D4] \n");
	//printf("\n");
	printf("TAP 0:    ");
	print256_i16(output[0]);
	printf("TAP 1:    ");
	print256_i16(output[1]);
	printf("TAP 2:    ");
	print256_i16(output[2]);
	printf("TAP 3:    ");
	print256_i16(output[3]);

#endif

}

//sort 4 channels of 4 pixels 4 times
void deinterleave4x4x4(__m256i* input, UINT16* output, int outstride)
{
#ifdef DEBUG
	printf("\n--------4 Way Interleaved DATA, 4 pixels, 4 times -----------------\n");
	printf("\n");
	printf("         [ ===========================================  YMM  =========================================   ] \n");
	printf("         [   TAP 1                   TAP 2                   TAP 3                   TAP 4               ] \n");
	printf("         [   A     B     C     D     A     B     C     D     A     B     C     D     A     B     C     D ] \n");

	for (int i = 0; i < 4; i++)
	{
		printf("Sample %d: ", i);
		print256_i16(input[i]);
	}

	printf("\n--------First Round insert into ZMM-----------------\n");
#endif

	//insert is faster than permute
	__m512i AB512 = _mm512_inserti32x8(_mm512_castsi256_si512(input[0]), input[1], 1);
	__m512i CD512 = _mm512_inserti32x8(_mm512_castsi256_si512(input[2]), input[3], 1);

#ifdef DEBUG
	print512_i16(AB512);
	print512_i16(CD512);

	printf("\n--------Resulting DATA in ZMM-----------------\n");
#endif

	__m512i idx = _mm512_set_epi16(61, 57, 53, 49, 45, 41, 37, 33, 29, 25, 21, 17, 13, 9, 5, 1, 60, 56, 52, 48, 44, 40, 36, 32, 28, 24, 20, 16, 12, 8, 4, 0);

	__m512i tmp = _mm512_mask_permutex2var_epi16(AB512, 0xFFFFFFFF, idx, CD512);



#ifdef DEBUG
	print512_i16(tmp);
#endif

	//save to 4 individual color channel frames spaced outstride apart in memory
	*(__m256i*)(output) = _mm512_extracti32x8_epi32(tmp, 0);
	*(__m256i*)(output + outstride) = _mm512_extracti32x8_epi32(tmp, 1);


	idx = _mm512_set_epi16(63, 59, 55, 51, 47, 43, 39, 35, 31, 27, 23, 19, 15, 11, 7, 3, 62, 58, 54, 50, 46, 42, 38, 34, 30, 26, 22, 18, 14, 10, 6, 2);

	tmp = _mm512_mask_permutex2var_epi16(AB512, 0xFFFFFFFF, idx, CD512);

#ifdef DEBUG
	print512_i16(tmp);
#endif

	*(__m256i*)(output + 2 * outstride) = _mm512_extracti32x8_epi32(tmp, 0);
	*(__m256i*)(output + 3 * outstride) = _mm512_extracti32x8_epi32(tmp, 1);




#ifdef DEBUG

	printf("\n");
	printf("         [ =============================================  YMM  ========================================= ] \n");
	printf("         [   Sample 1                Sample 2                Sample 3                Sample 4            ] \n");
	printf("         [   A     B     C     D     A     B     C     D     A     B     C     D     A     B     C     D ] \n");
	printf("TAP 0:    ");
	print256_i16(*(__m256i*)(output));
	printf("TAP 1:    ");
	print256_i16(*(__m256i*)(output + outstride));
	printf("TAP 2:    ");
	print256_i16(*(__m256i*)(output + 2 * outstride));
	printf("TAP 3:    ");
	print256_i16(*(__m256i*)(output + 3 * outstride));

#endif

}


//4 way interleaved data, 16 samples at a time (4 total permutes, 1 per output sample
//
void deinterleave16x4(__m256i* input, UINT16* output, int outstride)
{
#ifdef DEBUG
	printf("\n--------4 Way Interleaved DATA, 16 samples-----------------\n");

	for (int i = 0; i < 4; i++)
	{
		//    tosavearray[i] = _mm256_add_ph(tmp, _mm256_set1_ph((float)i));
		print256_i16(input[i]);
	}


	printf("\n--------First Round Merge-----------------\n");
#endif

	//take pairs of sequential values, skip two missing, then next two, ... from the YMM registers cast to ZMM
	__m512i idx10 = _mm512_set_epi16(47, 15, 46, 14, 45, 13, 44, 12, 43, 11, 42, 10, 41, 9, 40, 8, 39, 7, 38, 6, 37, 5, 36, 4, 35, 3, 34, 2, 33, 1, 32, 0);
	__m512i AB512 = _mm512_mask_permutex2var_epi16(_mm512_castsi256_si512(input[0]), 0xFFFFFFFF, idx10, _mm512_castsi256_si512(input[1]));
	__m512i CD512 = _mm512_mask_permutex2var_epi16(_mm512_castsi256_si512(input[2]), 0xFFFFFFFF, idx10, _mm512_castsi256_si512(input[3]));

#ifdef DEBUG
	print512_i16(AB512);
	print512_i16(CD512);

	printf("\n--------Resulting DATA-----------------\n");
	printf(" [===============================  Sample 0  ==================================================] [===============================  Sample 1  ==================================================]\n");
#endif

	//make the first half fully sequential by merging the above two vectors
	__m512i idx11 = _mm512_set_epi16(47, 46, 15, 14, 45, 44, 13, 12, 43, 42, 11, 10, 41, 40, 9, 8, 39, 38, 7, 6, 37, 36, 5, 4, 35, 34, 3, 2, 33, 32, 1, 0);
	//idx11 +16 to get the second half
	__m512i idx12 = _mm512_set_epi16(63, 62, 31, 30, 61, 60, 29, 28, 59, 58, 27, 26, 57, 56, 25, 24, 55, 54, 23, 22, 53, 52, 21, 20, 51, 50, 19, 18, 49, 48, 17, 16);
	__m512i AB512_2 = _mm512_mask_permutex2var_epi16(AB512, 0xFFFFFFFF, idx11, CD512);
	__m512i CD512_2 = _mm512_mask_permutex2var_epi16(AB512, 0xFFFFFFFF, idx12, CD512);

	//16-way deinterleaved
#ifdef DEBUG
	print512_i16(AB512_2);
	print512_i16(CD512_2);
#endif

	//save to 4 individual color channel frames spaced outstride apart in memory
	*(__m256i*)(output) = _mm512_extracti32x8_epi32(AB512_2, 0);
	*(__m256i*)(output + outstride) = _mm512_extracti32x8_epi32(AB512_2, 1);
	*(__m256i*)(output + 2*outstride) = _mm512_extracti32x8_epi32(CD512_2, 0);
	*(__m256i*)(output + 3*outstride) = _mm512_extracti32x8_epi32(CD512_2, 1);
}

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
