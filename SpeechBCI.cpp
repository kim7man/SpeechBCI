#include<stdio.h>
#include<stdlib.h>
#include<time.h>

#include "SpeechBCI.h"

#define N_Coeff 10
#define N_Data 37
#define N_Ch 30


void main()
{
	FILE *fp, *fout;
	int result;
	double sampleSeg[N_Ch][N_Data] ;
	bool eofFlag;


	fp = fopen("test_dataset.dat", "rt");
	fout = fopen("test_result.dat", "w");

	Classification test("model_R_vs_FN.dat", "model_F_vs_N.dat");
	eofFlag = false;
	do
	{
		int isExist;
		// read segment data until EOF occurs
		for(int seg_time=0; seg_time<N_Data; ++seg_time)
		{
			for(int ch_idx=0;ch_idx<N_Ch;++ch_idx)
			{
				double sampleInput;
				isExist = fscanf(fp, "%lf ", &sampleInput);
				sampleSeg[ch_idx][seg_time] = sampleInput;
				if(feof(fp))
				{
					eofFlag = true;
					break;
				}
			}
		}

		test.PutData((double*)sampleSeg, N_Ch, N_Data);
		result = test.ProcessData();
		printf("%d\n", result);
		fprintf(fout, "%d\n", result);
	}while(!eofFlag);

	fclose(fp);
	fclose(fout);
}
