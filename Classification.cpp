#include<memory.h>
#include<math.h>

#include<cstring>
using namespace std;

#include "Classification.h"


Classification::Classification(void)
{
	flagProcess = 0;

	//DLL initialize
	hFFTW = NULL;
}

Classification::Classification(char* modelName_base, char* modelName_fn)
{
	char tmp_base[30], tmp_fn[30];
	int l_base, l_fn;
	flagProcess = 0;

	headQue = 0;
	tailQue = 0;

	cntLabels = 0;
	labelCnt[0] = 0;
	labelCnt[1] = 0;
	labelCnt[2] = 0;
	iResult = 0;


	//load FFT dll file
	hFFTW = NULL;
	hFFTW = LoadLibrary(L"libfftw3x86.dll");
	if(hFFTW == NULL)
	{
		printf("0x%x\n", GetLastError());
	}

	// set function calls
	fftw_plan_r2r_1d = (pPlan1d)GetProcAddress(hFFTW, "fftw_plan_r2r_1d");
	fftw_execute = (pExecute)GetProcAddress(hFFTW, "fftw_execute");

	// initiallize FFT buffer
	planFFT = fftw_plan_r2r_1d(FFTSize, fftIn, fftOut, FFTW_R2HC, FFTW_MEASURE);
	idxSpect = 0;

	// initiallize spectrogram (ch x freq x time)
	fftSpect = new double**[NoChannel];
	for(int channel_index=0;channel_index<NoChannel;++channel_index)
	{
		fftSpect[channel_index] = new double*[(int)(FFTSize/2)];
		for(int frequency_index=0;frequency_index<(int)(FFTSize/2);++frequency_index)
		{
			fftSpect[channel_index][frequency_index] = new double[SpectSize];
		}
	}

	// load model from file
	model_base = svm_load_model(modelName_base);
	model_fn = svm_load_model(modelName_fn);

	strcpy(tmp_base, modelName_base);
	l_base = strlen(tmp_base)+2;
	strcpy(tmp_fn, modelName_fn);
	l_fn = strlen(tmp_fn)+2;
	for(int i=0;i<6; ++i)
	{
		tmp_base[l_base-i] = tmp_base[l_base-i-2];
		tmp_fn[l_fn-i] = tmp_fn[l_fn-i-2];
	}
	tmp_base[l_base-6] = '_';
	tmp_base[l_base-5] = 'c';
	tmp_fn[l_fn-6] = '_';
	tmp_fn[l_fn-5] = 'c';

	idx_base = svm_load_coeff(tmp_base);
	idx_fn = svm_load_coeff(tmp_fn);

	isFirst = true;
}

Classification::~Classification(void)
{
	delete model_base;
	delete model_fn;

	for(int i=0;i<NoChannel;++i)
	{
		for(int j=0;j<NoFreq;++j)
		{
			delete fftSpect[i][j];
		}
		delete fftSpect[i];
	}
	delete fftSpect;
}

bool Classification::isFinished()
{
	return flagProcess;
}

void Classification::PutData(double *inSegment, int nCh, int nCol)
{
	for(int channel_index=0; channel_index<nCh; ++channel_index)
	{
		if(headQue + nCol < QueueSize)
			memcpy(&queData[channel_index][headQue], &inSegment[channel_index*nCol], nCol*sizeof(double));
		else
		{
			int diff = QueueSize-headQue;
			memcpy(&queData[channel_index][headQue], &inSegment[channel_index*nCol], diff*sizeof(double));
			memcpy(&queData[channel_index][0], &inSegment[channel_index*nCol+diff], (nCol-diff)*sizeof(double));
		}
	}
	headQue += nCol;
	if(headQue >= QueueSize)
		headQue %= QueueSize;

	nChannel = nCh;
	flagProcess = 0;
};

int Classification::ProcessData()
{
	int diffHT;

	diffHT = (headQue + QueueSize - tailQue) % QueueSize;
	while(diffHT >= DataProcSize)
	{
		for(int channel_index=0; channel_index<nChannel; ++channel_index)
		{
			double inBuff[DataProcSize];

			// move moving window(data queue) to input buffer
			if(tailQue + DataProcSize <= QueueSize)
			{
				memcpy(inBuff, &queData[channel_index][tailQue], sizeof(double) * DataProcSize);
			}
			else
			{
				int nTail = QueueSize - tailQue;
				memcpy(inBuff, &queData[channel_index][tailQue], sizeof(double) * nTail);
				memcpy(&inBuff[nTail], &queData[channel_index][0], sizeof(double) * (DataProcSize-nTail));
			}

			// excute FFT
			memcpy(fftIn, inBuff, sizeof(double)* DataProcSize);
			fftw_execute(planFFT);

			// rearrange result FFT & append to spectrogram
			fftSpect[channel_index][0][idxSpect] = fftOut[0];
			for(int frequency_index=1;frequency_index<FFTSize/2;++frequency_index)
			{
				double Re = fftOut[frequency_index];
				double Im = fftOut[FFTSize-frequency_index];
				fftSpect[channel_index][frequency_index][idxSpect] = sqrt(Re*Re + Im*Im);
			}

		}

		++idxSpect;

		// when the data ready to calculate spectrogram
		if(idxSpect == SpectSize)
		{
			int cnt;
			int tmpResult;
			
		
			// the first time to make spectrogram -> from the start
			if(isFirst)
			{
				for(int channel_index=0;channel_index<nChannel;++channel_index)
				{
					for(int time_index=0;time_index<SpectSize;++time_index)
					{
						for(int frequency_index=0;frequency_index<NoFreq;++frequency_index)
						{
							double sum = 0;
							for(int freqBin_index=freqRange[frequency_index][0];freqBin_index<freqRange[frequency_index][1];++freqBin_index)
							{
								sum += fftSpect[channel_index][freqBin_index][time_index];
							}
							rdcSpect[channel_index][frequency_index][time_index] = sum / (freqRange[frequency_index][1] - freqRange[frequency_index][0] + 1);

							// accumulate spectral power along the time
							if(time_index == 0)
								accumSpect[channel_index][frequency_index][time_index] = rdcSpect[channel_index][frequency_index][time_index];
							else
								accumSpect[channel_index][frequency_index][time_index] = accumSpect[channel_index][frequency_index][time_index-1] + rdcSpect[channel_index][frequency_index][time_index];

							// 1000ms moving average 
							if(time_index >= TS_1000ms)
								accumSpect[channel_index][frequency_index][time_index] -= rdcSpect[channel_index][frequency_index][time_index-TS_1000ms];
						}
					}
					for(int frequency_index=0;frequency_index<NoFreq;++frequency_index)
					{
						double meanFeature = 0;
						for(int time_index=0;time_index<NoTime;++time_index)
						{
							featureMat[channel_index][frequency_index][time_index] = accumSpect[channel_index][frequency_index][TS_1000ms+(time_index*TS_150ms)] / TS_1000ms;
							meanFeature += featureMat[channel_index][frequency_index][time_index];
						}
						meanFeature /= NoTime;
						baseMat[channel_index][frequency_index] = meanFeature;

						for(int time_index=0;time_index<NoTime;++time_index)
						{
							normMat[channel_index][frequency_index][time_index] = log10(featureMat[channel_index][frequency_index][time_index] / meanFeature);
						}
					}
				}
				isFirst = false;
			}
			else	// second and after... -> calculate just for the shift
			{
				for(int channel_index=0;channel_index<nChannel;++channel_index)
				{
					for(int time_index=SpectSize-TS_150ms;time_index<SpectSize;++time_index)
					{
						for(int frequency_index=0;frequency_index<NoFreq;++frequency_index)
						{
							double sum = 0;
							for(int freqBin_index=freqRange[frequency_index][0];freqBin_index<freqRange[frequency_index][1];++freqBin_index)
							{
								sum += fftSpect[channel_index][freqBin_index][time_index];
							}
							rdcSpect[channel_index][frequency_index][time_index] = sum / (freqRange[frequency_index][1] - freqRange[frequency_index][0] + 1);

							// accumulate spectral power along the time
							accumSpect[channel_index][frequency_index][time_index] = accumSpect[channel_index][frequency_index][time_index-1] + rdcSpect[channel_index][frequency_index][time_index];
							// 1000ms moving average 
							accumSpect[channel_index][frequency_index][time_index] -= rdcSpect[channel_index][frequency_index][time_index-TS_1000ms];
						}
					}
					for(int frequency_index=0;frequency_index<NoFreq;++frequency_index)
					{
						featureMat[channel_index][frequency_index][NoTime-1] = accumSpect[channel_index][frequency_index][TS_1000ms+((NoTime-1)*TS_150ms)] / TS_1000ms;
						normMat[channel_index][frequency_index][NoTime-1] = log10(featureMat[channel_index][frequency_index][NoTime-1] / baseMat[channel_index][frequency_index]);
					}
				}				
			}


			// linearize features (3d->1d, Ch x Time x Freq)
			nFeatures_base = nChannel * NoTime * NoFreq;
//			nFeatures_base = idx_base->l;
//			nFeatures_fn = idx_fn->l;
//			featureMat_base = new struct svm_node [nFeatures_base + 1];
//			featureMat_fn = new struct svm_node [nFeatures_fn + 1];
			featureMat_base = new struct svm_node [nFeatures_base + 1];

			cnt = 0;
			for(int channel_index=0; channel_index<nChannel; ++channel_index)
			{
				for(int freq_index=0; freq_index<NoFreq; ++freq_index)
				{
					for(int time_index=0; time_index<NoTime; ++time_index)
					{
						featureMat_base[cnt].index = cnt+1;
		//				featureMat_base[cnt++].value = normMat[idx_base->coef[idx][0]][idx_base->coef[idx][1]][idx_base->coef[idx][2]];
						featureMat_base[cnt++].value = featureMat[channel_index][freq_index][time_index];
					}
				}
			}
			featureMat_base[nFeatures_base].value = -1;
/*
			cnt = 0;
			for(int idx=0; idx<nFeatures_base; ++idx)
			{
				featureMat_base[cnt].index = cnt+1;
//				featureMat_base[cnt++].value = normMat[idx_base->coef[idx][0]][idx_base->coef[idx][1]][idx_base->coef[idx][2]];
				featureMat_base[cnt++].value = featureMat[idx_base->coef[idx][0]][idx_base->coef[idx][1]][idx_base->coef[idx][2]];
			}
			featureMat_base[nFeatures_base].value = -1;

			cnt = 0;
			for(int idx=0; idx<nFeatures_fn; ++idx)
			{
				featureMat_fn[cnt].index = cnt+1;
//				featureMat_fn[cnt++].value = normMat[idx_fn->coef[idx][0]][idx_fn->coef[idx][1]][idx_fn->coef[idx][2]];
				featureMat_fn[cnt++].value = featureMat[idx_fn->coef[idx][0]][idx_fn->coef[idx][1]][idx_fn->coef[idx][2]];
			}
			featureMat_fn[nFeatures_fn].value = -1;
*/
			// classification
			// test if the block is baseline or not
			tmpResult = (int)svm_predict(model_base, featureMat_base);
			if(tmpResult)
			{
				// differentiate face/number
				tmpResult = (int)svm_predict(model_fn, featureMat_base);
//				tmpResult = (int)svm_predict(model_fn, featureMat_fn);
			}

			++labelCnt[tmpResult];
			if(cntLabels < NoLabel-1)
			{// voting initiallize
				recentLabels[cntLabels++] = tmpResult;
				iResult = 0;
			}
			else
			{
				int labelMax = 0;

				recentLabels[cntLabels] = tmpResult;
				
				// voting system (more than half, default : neutral) 
				if(labelCnt[labelMax] < labelCnt[1])
					labelMax = 1;
				if(labelCnt[labelMax] < labelCnt[2])
					labelMax = 2;
				if(labelCnt[labelMax] < NoLabel/2)
					labelMax = 0;
				iResult = labelMax;

				// shift labels for next use
				labelCnt[recentLabels[0]]--;
				for(int label_index=0;label_index<NoLabel-1;++label_index)
				{
					recentLabels[label_index] = recentLabels[label_index+1];
				}
			}

			delete featureMat_base;
//			delete featureMat_fn;
			
			
			// shift 150ms
			for(int channel_index=0;channel_index<nChannel;++channel_index)
			{
				for(int frequency_index=0;frequency_index<(int)(FFTSize/2);++frequency_index)
				{
					for(int time_index=0;time_index<SpectSize-TS_150ms;++time_index)
					{
						fftSpect[channel_index][frequency_index][time_index] = fftSpect[channel_index][frequency_index][time_index+TS_150ms];
					}
				}
				for(int frequency_index=0;frequency_index<NoFreq;++frequency_index)
				{
					for(int time_index=0;time_index<SpectSize-TS_150ms;++time_index)
					{
						rdcSpect[channel_index][frequency_index][time_index] = rdcSpect[channel_index][frequency_index][time_index+TS_150ms];
						accumSpect[channel_index][frequency_index][time_index] = accumSpect[channel_index][frequency_index][time_index+TS_150ms];
					}
					for(int time_index=0;time_index<NoTime-1;++time_index)
					{
						featureMat[channel_index][frequency_index][time_index] = featureMat[channel_index][frequency_index][time_index+1];
					}
				}
			}
			idxSpect -= TS_150ms;
		}

		// shift moving window
		tailQue += T_50ms;
		if(tailQue >= QueueSize)
			tailQue %= QueueSize;
		diffHT = (headQue + QueueSize - tailQue) % QueueSize;
	}


	CString a;

	flagProcess = 1;

	return iResult;
}



