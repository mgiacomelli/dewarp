#include "stdafx.h"
#include "stdio.h"
#include <wtypes.h>
#include <math.h>
#include <malloc.h >

#include <stdlib.h>
#include <windows.h>

#include <immintrin.h>
#include "avx_print.h"
#include "vec_permute.h"

#define PI 3.14159265358979323846 

void saveBufferToDisk(char* buf, int length);
void saveBufferToDisk(char* buffer, int length, int num);
void wsfiltgen(float* fir, int nt, double bw);
void fir_avx512_fp16(int taps_samp, __bfloat16* barrayh, __m256i* startPtr, __m256h* filtered_data, int num_bits, int dc_offset);
void fir_avx256_fp16(int taps_samp, __bfloat16* barrayh, __m256i* startPtr, __m256h* filtered_data, int num_bits, int dc_offset);
void cubicHermitePolyFp16(__m256i* input, __m256i* output, __m256h frac_scalar_h, bool renormalize);

//VS2022 is missing this intrinsic for some reason
#if 1 &&  _MSC_VER && !__INTEL_COMPILER	
__m256h _mm256_set1_ph(float f)
{
	return _mm512_cvtps_ph(_mm512_set1_ps(f), 0);
}
#endif

//the DC offset for the signal
const short dc_offset = 32000;
//number of bits in the ADC
const short num_bits = 14;

int main()
{

	const int channels = 16;
	const int sampleSize = 32;	//channels times bytes per sample

	//resonant scanner parameters
	const int samples = 6656; //bidirectional, so really 2x3328
	const int sampleClockLength = 6700;	//the number of ADC samples between triggers (Sampling_rate/resonant_scanner_frequency)
	//note that there are usually a few samples that are not recorded between triggers to allow for the resonant scanner clock drift

	//the alignment between the forward and backwards sweeps
	const float shift = 34;

	//check that there are not more samples per buffer than the duration of one bidirectional sweep
	if (sampleClockLength < samples)
	{
		printf("ERROR:  number of samples per buffer (%d) must be less than or equal to the scanner period in samples (%d)!\n", samples, sampleClockLength);
		return -1;
	}

	const int galvoResetLines = 0;	//any non-image lines while the galvo moves into position to begin the frame
	const int lines = 1024 + galvoResetLines;
	const int outputsamples = 2048;	//how many samples to reconstruct per line after interpolation

	const bool usefp16 = true;		//required for this version
	const bool useAVX512 = true;	//required for this version

	//Do FIR filtering inline (faster) or before hand (slower)
	const bool useFirPreFiltering = false;
	const bool useFirInlineFiltering = true;

	//to improve benchmark accuracy
	const int timesToAverage = 40;

	//allocate a buffer worth of input data, force 64 bit alignment for ZMM registers, add 4096 bytes before and after for the FIR filter
	UINT16* data_p = (UINT16*)_aligned_malloc(samples * lines * sampleSize + 4096*2, 64); //6656 samples per trigger, 528 triggers per frame, 2 16 bit channels

	//aligned malloc does not initialize to zero and there is no aligned calloc, and should initialize accounting for DC offset
	UINT16 zeroVal = dc_offset;

	for (int i = 0; i < 2048; i++)
	{
		data_p[i] = zeroVal;
		data_p[samples * lines * channels + 2048 + i] = zeroVal;
	}

	//hack so that we can safely access the addresses corresponding to missing samples at the start (and end) if no flyback samples
	UINT16* data = data_p + 2048;

	//allocate a frame of output data
	UINT16* output = (UINT16*)_aligned_malloc(outputsamples * (lines - galvoResetLines) * 2 * sizeof(UINT16) * channels, 64); //6656 samples per trigger, 528 triggers per frame, 2 16 bit channels


	//get a buffer worth of ADC data
	FILE* pFile;
	char buf[255];
	if (sampleSize == 4)
	{	
		//not implemented in this version of the code

	}
	else if (sampleSize == 3)
	{

	}
	else if (sampleSize == 32)
	{

		sprintf(buf, "warpedFrame_%dchan_2048.bin", 16);	//2048x2048 test 
		pFile = fopen(buf, "rb");
		if (pFile == NULL)
		{
			printf("Could not open file: %s\n", buf);
			return -2;
		}

		int bytesRead = fread(data, 1, samples * lines * sampleSize, pFile);
		printf("read %d bytes\n", bytesRead);

		if (bytesRead != samples * lines * sampleSize)
		{
			perror("Error: ");

			return -1;
		}
	}

	fclose(pFile);

	//calculate the dewarping indecies

	float index_table[outputsamples * 2];

	//match variable names from the paper
	int m = sampleClockLength;
	int n = outputsamples;
	for (int i = 0; i < n; i++)
	{
		index_table[i] = (float) m / 2.0 / PI * acos(1 - 2.0 * (i + 0.5) / n) - shift;

		//make sure that all 4 filter taps (3 to the left of the index) will be in bounds
		if (index_table[i] < 3)
			index_table[i] = 3;
	}
	//backwards sweep
	for (int i = n; i < 2 * n; i++)
	{
		index_table[3 * n - i - 1] = (float) m / 2.0 + m / 2.0 / PI * acos(2.0 / n * (i + 0.5 - n) - 1) - shift;

		//make sure we don't use samples after the end of the data 
		if (index_table[3 * n - i - 1] > samples - 1)
			index_table[3 * n - i - 1] = samples - 1;

		//HACK so that we don't break our in place FIR filtering, TODO:  fix this 
		if (index_table[3 * n - i - 1] < sampleClockLength / 2 + 3)
			index_table[3 * n - i - 1] = sampleClockLength / 2 + 3;
	}

	saveBufferToDisk((char*)index_table, outputsamples * 2 * sizeof(float), 10);

	/*begin on the FIR calculation section*/

	//%calculate the cutoff frequency assuming that the image is Nyquist sampled at index samples / 4 (center FOV)
	double oversample[outputsamples];
	for (int i = 0; i < outputsamples - 1; i++)
		oversample[i] = index_table[i + 1] - index_table[i];

	//first index and end are missing / wrong due to diff, so duplicate 2 and end - 1
	oversample[0] = oversample[1];
	oversample[outputsamples - 1] = oversample[outputsamples - 2];

	//cutoff frequencies, note symetric so can use for both directions
	double cutoff[outputsamples];
	for (int i = 0; i < outputsamples; i++)
		cutoff[i] = 1 / oversample[i];

	//calculate the number of taps for each pixel assuming the transition band will be 3.5 / Ntaps wide
	int outtaps[outputsamples];
	double shifted_cutoff[outputsamples];
	for (int i = 0; i < outputsamples; i++)
	{
		//initial guess
		int ntaps = (int) 3 * ceil(oversample[i]);

		//calculate the transistion band half width
		double tband = 1.75 / ntaps;

		// make sure the tband isn't too large to be useful
		int newtaps = ntaps;
		while (tband > 0.25)
		{
			newtaps = newtaps + 1;
			tband = 1.75 / newtaps;
		}

		// first shift the cut off by half the transition band
		shifted_cutoff[i] = cutoff[i] + tband;

		//clamp cutoff at 1
		if (shifted_cutoff[i] > 1)
			shifted_cutoff[i] = 1;

		//always use even n for integer group delay(fir1 adds 1)
		if (newtaps % 2 == 1)
			newtaps = newtaps + 1;
		outtaps[i] = newtaps;
	}

	//we have taps and cutoffs in terms of output pixels, but will be filtering the raw data in samples
	//nearest neighbor the taps/cutoff to the sample domain from pixel domain
	float ii = index_table[0];
	int tablepos = 0;


	float cutoff_samp[sampleClockLength / 2];
	int taps_samp[sampleClockLength / 2];

	//get any samples on the edge
	for (int i = 0; i < ii; i++)
	{
		cutoff_samp[i] = shifted_cutoff[0];
		taps_samp[i] = outtaps[0];
	}

	for (int i = (int) ii; i < sampleClockLength / 2 + 1; i++)	//we got to <= and address i-1 to match the matlab code, might not be ideal
	{
		if (i >= ii)
		{
			tablepos++;
			ii = index_table[tablepos];
		}

		if (tablepos >= outputsamples)
		{
			//points after the last table value aren't used anyway
			cutoff_samp[i - 1] = shifted_cutoff[outputsamples - 1];
			taps_samp[i - 1] = outtaps[outputsamples - 1];
		}
		else
		{
			cutoff_samp[i - 1] = shifted_cutoff[tablepos];
			taps_samp[i - 1] = outtaps[tablepos];
		}

	}

	saveBufferToDisk((char*)shifted_cutoff, outputsamples * sizeof(double), 103);

	saveBufferToDisk((char*)cutoff_samp, sampleClockLength / 2 * sizeof(float), 501);

	int barrayMaxCoefficients = 100000;
	float* barray = (float*)malloc(barrayMaxCoefficients * sizeof(float));

	int bpos = 0;
	// calculate the filter coefficients
	for (int i = 0; i < samples / 2; i++)	//note that sampleClockLength >= samples so not all of these are used
	{
		if (cutoff_samp[i] > 0.99)
		{
			taps_samp[i] = 0;	//rounds up to 1
			barray[bpos++] = 1;

		}
		else
		{

			wsfiltgen(&barray[bpos], taps_samp[i] + 2, cutoff_samp[i]);	//N+2 to match matlab fir1, there are actually N+1 total taps returned

			//it may actually be better to use old data since otherwise there is a change in group delay for the first ~20 samples or so
			if (0)
			{
				if (i < taps_samp[i] / 2)
				{
					//handle the boundary condition
					int tozero = taps_samp[i] / 2 - i;
					for (int j = 0; j < tozero; j++)
					{
						barray[bpos + j] = 0;
					}

					//renormalize so DC is correct
					float sum = 0;
					for (int j = 0; j < taps_samp[i] + 1; j++)
					{
						sum += barray[bpos + j];
					}
					for (int j = 0; j < taps_samp[i] + 1; j++)
					{
						barray[bpos + j] /= sum;
					}
				}
			}

			bpos += (taps_samp[i] + 1);	//returns N+1 taps

		}

	}
	printf("Using %d total FIR taps per line/channel (%d per frame)\n", bpos + 1, (bpos + 1) * (lines - galvoResetLines) * 2 * channels);
	__bfloat16* barrayh = (__bfloat16*)malloc((bpos + 9) * sizeof(__bfloat16));

	for (int i = 0; i < bpos + 1; i += 8)
	{
		//VS2022 scalar intrinsics for FP16 are missing, so do this with intrinsics so it can compile with VS or Clang
		__m256 b = _mm256_loadu_ps(&(barray[i]));

		__m128h bh = _mm256_cvtxps_ph(b);

		_mm_store_ph(&(barrayh[i]), bh);

	}

	saveBufferToDisk((char*)taps_samp, sampleClockLength / 2 * sizeof(int), 102);
	saveBufferToDisk((char*)barray, bpos * sizeof(float), 100);
	saveBufferToDisk((char*)barrayh, bpos * sizeof(__bfloat16), 99);


	//unpack and dewarp the ADC data into pixels

	int startfpos = samples * sampleSize * galvoResetLines;
	int outsamps = 0;

	//do not assume a square dataset, can be 2048x1024 for example 
	const int outstride = outputsamples * 2 * (lines - galvoResetLines);	//width of one spectral channel worth of output in pixels, 2 is because each sweep gives a forward and a backwards line

	//how many samples we will process in one batch (16 32 byte samples, resulting in 256 bit writes to each of 16 output channel buffers)
	const int sampleStepSize = 16;

	//check for invalid output sample size
	if (outputsamples % sampleStepSize != 0)
	{
		printf("Output samples (%d) needs to be divisible by sampleStepSize (%d)!\n", outputsamples, sampleStepSize);
		return 0;
	}

	//buffer for FIR filtering, hold only one direction of a sweep at a time
	__m256* filtered_data = (__m256*)_aligned_malloc(2 * sampleClockLength / 2 * sampleSize, 64);	//2 because we go from int16 to float samples32)
	UINT8* filtered_data_ptr;

	if (filtered_data == NULL)
		return -1;

	LARGE_INTEGER frequency;
	LARGE_INTEGER start;
	LARGE_INTEGER end;
	double interval;

	QueryPerformanceFrequency(&frequency);
	QueryPerformanceCounter(&start);


	for (int averages = 0; averages < timesToAverage; averages++)	//only used to average for better benchmarks
	{
		outsamps = 0;
		startfpos = samples * sampleSize * galvoResetLines;

		for (int i = galvoResetLines; i < lines; i++)			//slow axis loop
		{
			//forward line 
			UINT8* linedata = ((UINT8*)data) + startfpos;
			filtered_data_ptr = linedata;

			//FIR filtering loop
			if (useFirPreFiltering)
			{
				int count = 0;
				
				for (int ii = 0; ii < samples / 2; ii++)	//I think we don't need to filter past samples/2 since the interpolation will never access points beyond this
				{

					__m128i* startPtr = (__m128i*) (linedata + (ii - taps_samp[ii] / 2) * sampleSize);
					//if(ii==0)printf("%p %p %d currentTaps: %d\n", linedata, startPtr, (ii - taps_samp[ii] / 2)* sampleSize, taps_samp[ii]);
					fir_avx512_fp16(taps_samp[ii], barrayh + count, (__m256i*)startPtr, (__m256h*)filtered_data + ii, num_bits, dc_offset);
					count += taps_samp[ii] + 1;
				}

				filtered_data_ptr = (UINT8*)filtered_data;

				if (0 && i == 0) {
					saveBufferToDisk((char*)filtered_data, samples / 2 * sampleSize * 2, 400);	//FP32 or FP16
					saveBufferToDisk((char*)linedata, sampleClockLength / 2 * sampleSize, 401);		//int16, figure;plot(d(1:16:end)) to see a fast line
				}
				//note:  output is very slightly different (+/- 1) vs matlab.  Can test with:
				/*f=fopen('myDebugData400.bin');
				d=fread(f, 'float');
				fclose(f); figure;plot(d(1:16:end)); hold on; ;plot(filtered(1:dataLen/2))*/

			}
			else if (useFirInlineFiltering)
			{
				filtered_data_ptr = (UINT8*)filtered_data;
			}
			else
			{
				//no filtering
				filtered_data_ptr = linedata;
			}
			int currentFilterPos = 0;
			int count = 0;
			for (int j = 0; j < outputsamples; j += sampleStepSize)
			{

				float frac;
				__m256i tosavearray[sampleStepSize];

				//outer step size loop
				for (int kk = 0; kk < sampleStepSize; kk += 4)
				{
					//inner 4 element step size loop
					for (int k = 0; k < sampleStepSize / 4; k++)
					{

						int pos = floor(index_table[j + k + kk]);
						frac = index_table[j + k + kk] - floor(index_table[j + k + kk]);

						//printf("pos/frac: %d %f\n", pos, frac);
						__m256h frac_scalar_h = _mm256_set1_ph(frac);

						__m256i data[4];

						if (usefp16)
						{
							if (useAVX512)
							{
								int tmp = (pos)*sampleSize;


								//it doesn't make sense to call it on zero, one... directly since we'll compute many duplicates
								//instead, try filling filtered_data here as we get to points for locality and skipping of unused points

								if (useFirInlineFiltering) {

									//skip unused points
									while (currentFilterPos < pos - 3)
									{
										count += (taps_samp[currentFilterPos] + 1);
										currentFilterPos++;
									}

									//the compiled assembly here looks ugly, see if it can be cleaned up by simplifying pointers
									for (; currentFilterPos <= pos; currentFilterPos++)
									{
										int currentTaps = taps_samp[currentFilterPos];
										__m256i* startPtr = (__m256i*) (linedata + (currentFilterPos - currentTaps / 2) * sampleSize);
										fir_avx512_fp16(currentTaps, &barrayh[count], startPtr, (__m256h*)&(filtered_data_ptr[currentFilterPos * sampleSize]), num_bits, dc_offset);
										count += currentTaps + 1;
									}
								}

								//for now, load just 16 samples (256 bits/32 bytes) rather than try packing two samples into ZMM registers
								//TODO:  try doing 512 bits per pass
								data[0] = _mm256_load_ph((float*)&(filtered_data_ptr[tmp]));
								tmp -= sampleSize;

								data[1] = _mm256_load_ph((float*)&(filtered_data_ptr[tmp]));
								tmp -= sampleSize;

								data[2] = _mm256_load_ph((float*)&(filtered_data_ptr[tmp]));
								tmp -= sampleSize;

								data[3] = _mm256_load_ph((float*)&(filtered_data_ptr[tmp]));

							}
							else
							{
								//AVX2 not implemented in this version
							}
						}
						else
						{
							//load from filtered_data which is already float32 using _mm256_load_ps()
							//removed FP32 from this version

						}

						cubicHermitePolyFp16(data, &tosavearray[k + kk], frac_scalar_h, !useFirInlineFiltering && !useFirPreFiltering);

						//finished 1 16 channel pixel

					}

					//finished 4 16 channel pixels
					//in fp16, now contains 4 channel 16 bit integers arranged in YMM [ABCDEFGHIJKLMNOP][ABCDEFGHIJKLMNOP]...

				}

				//finished 16, 16 channel pixels, now deinterleave them into individual buffers outstride apart in memory
				deinterleave16x16(tosavearray, (output + outsamps), outstride);

				outsamps += sampleStepSize;

			}
			//printf("outsamps: %d", outsamps);
			startfpos += samples * sampleSize;

			//backwards line

			//set the currentFilterPos and count to be the last samples of the buffers since we walk the reverse scan backwards 
			currentFilterPos = (samples / 2 - 1);
			//currentFilterPos = (samples- 1);
			//need index of hte start of the last set of filter taps
			count = bpos;

			//FIR filtering loop #2
			if (useFirPreFiltering)
			{
				int count = 0;
				//fow now just filter the entire line, eventually should skip the samples we never use (not too many of these at 2048 but more significant at 1024)
				//for (i = 0; i < sampleClockLength / 2; i++)	//note that we filter to sampleClockLength/2 (6700/2), which is the turnaround point, the backscan will be from here to 'samples' (6656)
				for (int ii = 0; ii < samples / 2; ii++)	//I think we don't need to filter past samples/2 since the interpolation will never access points beyond thisuyh
				{
					//int startindex = i - taps_samp[i] / 2;
					__m128i* startPtr = (__m128i*) (linedata + (sampleClockLength / 2 + ii - taps_samp[ii] / 2) * sampleSize);	//shift by samples / 2 because this is the second half of the scan
					//could fold this into the previous step, should benchmark and see if its faster to do both halve together, might increase cache misses though
					fir_avx512_fp16(taps_samp[ii], barrayh + count, (__m256i*)startPtr, (__m256h*)filtered_data + ii, num_bits, dc_offset);
					//if (ii == 0)printf("%p %p %d currentTaps: %d\n", linedata, startPtr, (sampleClockLength / 2 + ii - taps_samp[ii] / 2) * sampleSize, taps_samp[ii]);
					/*if (ii == samples / 2 - 1 - 1)
					{
						printf("fir_avx512_fp16(\n %d\n %p\n %p\n %p\n)\n", taps_samp[ii], barrayh + count, startPtr, filtered_data + ii);
						printf("count: %d\n", count);
					}*/
					count += taps_samp[ii] + 1;
				}
				filtered_data_ptr = (UINT8*)filtered_data;

				if (0) {
					saveBufferToDisk((char*)filtered_data, samples / 2 * sampleSize * 2, 402);	//floats
					//saveBufferToDisk((char*)linedata, sampleClockLength / 2 * sampleSize, 403);		//int16, figure;plot(d(1:16:end)) to see a fast line
					saveBufferToDisk((char*)linedata, samples * sampleSize, 403);		//figure; plot(d2(1:16 : end)); hold on; plot(6656 / 2 + [1:length(d) / 16], d(1:16 : end))
				}
			}
			else if (useFirInlineFiltering)
			{
				filtered_data_ptr = (UINT8*)filtered_data;
			}
			else
			{
				//no filtering, just point to the position in the raw data
				filtered_data_ptr = (linedata + (sampleClockLength / 2) * sampleSize);;
				//filtered_data_ptr = linedata;
			}

			for (int j = outputsamples - 1; j >= 0; j -= sampleStepSize)
			{

				float frac;
				__m256i tosavearray[sampleStepSize];

				//outer step size loop
				for (int kk = 0; kk < sampleStepSize; kk += 4)
				{
					//inner 4 element step size loop
					for (int k = 0; k < sampleStepSize / 4; k++)
					{

						int pos = (int)floor(index_table[outputsamples + j - k - kk]) - sampleClockLength / 2;
						frac = index_table[outputsamples + j - k - kk] - (float)floor(index_table[outputsamples + j - k - kk]);

						__m256h frac_scalar_h = _mm256_set1_ph(frac);

						__m256i data[4];


						if (usefp16)
						{
							if (useAVX512)
							{
								//load from filtered_data which is in fp16

								//512 bits (64 bytes) is two whole samples in fp16
								//shift tmp over samples / 2  since we're only storing the second half of the filtered data now (delete if we compute it all at once in the future)
								int tmp = (pos)*sampleSize;	//can probably simplify this once in situ works

								//printf("tmp: %d currentFilterPos: %d\n", tmp, currentFilterPos);

								if (useFirInlineFiltering) {
									//skip any samples we won't use
									while (currentFilterPos > pos)
									{
										count -= (taps_samp[currentFilterPos] + 1);
										currentFilterPos--;
										//printf("skipped (%d) because pos: %d index_table: %f\n", currentFilterPos, pos, index_table[outputsamples + j - k - kk]);
									}

									//walk memory backwards
									for (; currentFilterPos >= pos - 3; currentFilterPos--)
									{
										int currentTaps = taps_samp[currentFilterPos];
										count -= (currentTaps + 1);


										//printf("currentFilterPos: %d index_table: %f currentTaps: %d\n", currentFilterPos, index_table[outputsamples + j - k - kk], currentTaps);
										//todo:  optimize by precomputing
										__m256i* startPtr = (__m256i*) (linedata + (currentFilterPos - currentTaps / 2 + sampleClockLength / 2) * sampleSize);
										//if(currentFilterPos==0)printf("%p %p offset: %d currentTaps: %d pos: %d\n", linedata, startPtr, (currentFilterPos / 32 - currentTaps / 2)* sampleSize, currentTaps, tmp);
										/*if (currentFilterPos == sampleSize * (samples / 2 - 1 - 1))
										{
											printf("!fir_avx512_fp16(\n %d\n %p\n %p\n %p\n)\n", currentTaps, &barrayh[count], startPtr, &(filtered_data_ptr[currentFilterPos]));
											printf("count: %d\n", count);
											//exit(0);
										}*/

										fir_avx512_fp16(currentTaps, &barrayh[count], startPtr, (__m256h*)&(filtered_data_ptr[currentFilterPos * sampleSize]), num_bits, dc_offset);

									}
								}
								

								data[0] = _mm256_load_ph((float*)&(filtered_data_ptr[tmp]));
								tmp -= sampleSize;

								data[1] = _mm256_load_ph((float*)&(filtered_data_ptr[tmp]));
								tmp -= sampleSize;

								data[2] = _mm256_load_ph((float*)&(filtered_data_ptr[tmp]));
								tmp -= sampleSize;

								data[3] = _mm256_load_ph((float*)&(filtered_data_ptr[tmp]));


							}
							else
							{
								//load from filtered_data which is in fp16

							}

						}

						else
						{

						}
						cubicHermitePolyFp16(data, &tosavearray[k + kk], frac_scalar_h, !useFirInlineFiltering && !useFirPreFiltering);

					}

					//finished sampleStepSize (4) 16 channel pixels

					//tosavearray now contains 32 bit integers arranged in i256 [ABCDEFGH_0][IJKLMNOP_0][ABCDEFGH_1][IJKLMNOP_1]..., convert to 16 bit and then write out in pairs

				}


				//finished 16, 16 channel pixels, now deinterleave them into individual buffers outstride apart in memory
				deinterleave16x16(tosavearray, (output + outsamps), outstride);

				outsamps += sampleStepSize;
			}

		}
	}

	QueryPerformanceCounter(&end);
	
	//save the processed frame to disk
	saveBufferToDisk((char*)output, outsamps * sizeof(UINT16) * channels, 0);

	interval = (double)(end.QuadPart - start.QuadPart) / frequency.QuadPart;

	outsamps *= timesToAverage;

	printf("%f (%f fps)\n", interval, timesToAverage / interval);

	printf("processed %d samples of %d channels (%d MP)\n", outsamps, channels, outsamps / 1024 / 1024 / timesToAverage);
	printf("Total throughput:  %f Mspectra per second\n", (outsamps / 1024.0 / 1024.0) * channels / interval);

	//free buffers
	_aligned_free(data_p);
	_aligned_free(output);
	_aligned_free(filtered_data);

	return 0;
}

void cubicHermitePolyFp16(__m256i* input, __m256i* output, __m256h frac_scalar_h, bool renormalize)
{

	__m256h zero = input[0];
	__m256h one = input[1];
	__m256h two = input[2];
	__m256h three = input[3];

	if (renormalize)
	{
		//only used if no FIR, the FIR filter handles normalization otherwise

		//shift off any zeros in LSBs
		zero = _mm256_srli_epi16(zero, (short)(16 - num_bits));
		one = _mm256_srli_epi16(one, (short)(16 - num_bits));
		two = _mm256_srli_epi16(two, (short)(16 - num_bits));
		three = _mm256_srli_epi16(three, (short)(16 - num_bits));

		//remove any dc offset, critical for fp16 due to limited dynamic range
		zero = _mm256_sub_epi16(zero, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));
		one = _mm256_sub_epi16(one, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));
		two = _mm256_sub_epi16(two, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));
		three = _mm256_sub_epi16(three, _mm256_set1_epi16(dc_offset >> (16 - num_bits)));

		//convert to fp16 (signed)
		zero = _mm256_cvtepi16_ph(zero);
		one = _mm256_cvtepi16_ph(one);
		two = _mm256_cvtepi16_ph(two);
		three = _mm256_cvtepi16_ph(three);

	}

	//float tmp0 = (x2 - x1);				
	__m256h tmp = _mm256_sub_ph(two, one);

	//float dx = x0 - x3;
	__m256h diff = _mm256_sub_ph(zero, three);

	//acc0 = (3 * tmp0 + dx) / 2;
	__m256h acc = _mm256_fmadd_ph(_mm256_set1_ph(3.0f), tmp, diff);
	acc = _mm256_mul_ph(acc, _mm256_set1_ph(0.5f));

	//acc0 = (acc0 * frac);
	acc = _mm256_mul_ph(acc, frac_scalar_h);

	//tmp0 = (5 * x2 + x0);
	tmp = _mm256_fmadd_ph(_mm256_set1_ph(5.0f), two, zero);

	//float twoSumx = 2 * x1 + x3;
	__m256h twoSum = _mm256_fmadd_ph(_mm256_set1_ph(2.0f), one, three);

	//acc0 = acc0 + twoSumx - (tmp0 / 2);
	tmp = _mm256_fmadd_ph(_mm256_set1_ph(-0.5f), tmp, twoSum);
	acc = _mm256_add_ph(acc, tmp);

	//acc0 = (acc0 * frac);
	acc = _mm256_mul_ph(acc, frac_scalar_h);

	//acc0 = acc0 + (x1 - x3) / 2;
	diff = _mm256_sub_ph(one, three);
	acc = _mm256_fmadd_ph(diff, _mm256_set1_ph(0.5f), acc);

	//acc0 = (acc0 * frac);
	acc = _mm256_mul_ph(acc, frac_scalar_h);

	//acc0 = acc0 + x2;
	acc = _mm256_add_ph(acc, two);

	//in 256h we now have 16 complete channels processed for 1 sample

	*output = _mm256_cvtph_epi16(acc);

}


void saveBufferToDisk(char* buffer, int length)
{
	saveBufferToDisk(buffer, length, 0);
}

void saveBufferToDisk(char* buffer, int length, int num)
{
	FILE* pFile;
	char buf[255];
	sprintf(buf, "myDebugData%d.bin", num);
	pFile = fopen(buf, "wb");
	fwrite(buffer, sizeof(char), length, pFile);
	fclose(pFile);

}

void wsfiltgen(float* fir, int nt, double bw)
{
	double sum = 0;
	for (int i = 1; i <= nt - 1; i++)
	{
		double  a = (i - nt / 2) * PI * bw;
		double ys;
		if (a == 0)
			ys = 1;
		else
			ys = sin(a) / a;

		double yw = 0.54 - 0.46 * cos((i - 1) * 2.0 * PI / (nt - 2));

		fir[i - 1] = (float)yw * bw * ys;
		//printf("%f ", fir[i-1]);
		sum += fir[i - 1];
	}

	//normalize intensity
	for (int i = 0; i < nt - 1; i++)
		fir[i] /= sum;
}