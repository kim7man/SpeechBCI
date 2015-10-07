#include<memory.h>
#include<math.h>
#include<string.h>

#include<cstring>
using namespace std;

#include "Classification.h"


const double a[13] = {1,-6.74986612434222, 20.6816848333804, -38.6350580503587, 49.8950151033794, -47.5378229326030, 34.3183052306746, -18.7995025441968, 7.73677111762930, -2.34340252018987, 0.492840603242827, -0.0630345108008912, 0.00406980731467992};
const double b[13] = {0.0121493449463388, 0, -0.0728960696780330, 0, 0.182240174195082, 0, -0.242986898926777, 0, 0.182240174195082, 0, -0.0728960696780330, 0, 0.0121493449463388};


Classification::Classification(void)
{
	flagProcess = 0;

	//DLL initialize
	hFFTW = NULL;
}

Classification::Classification(char* modelName_base, char* modelName_fn)
{
	char filename_base[50], filename_fn[50];
	char *tmpfilename_base, *tmpfilename_fn;
	int l_base, l_fn;
	FILE *fp_base, *fp_fn;
	flagProcess = 0;

	headQue = 0;
	tailQue = 0;

	cntLabels = 0;
	labelCnt[0] = 0;
	labelCnt[1] = 0;
	labelCnt[2] = 0;
	iResult = 0;

	// initialize freq bin
	for(int index_freq=0; index_freq<12; ++index_freq)
	{
		freqRange[index_freq][0] = index_freq*4;
		freqRange[index_freq][1] = (index_freq+1)*4-1;
	}
	//for(int index_freq=0; index_freq<49; ++index_freq)
	//{
	//	freqRange[index_freq][0] = index_freq;
	//	freqRange[index_freq][1] = index_freq;
	//}

	//set filter
	filter = new Filtfilt*[NoChannel];
	for(int idxChannel=0; idxChannel<NoChannel; ++idxChannel)
	{
		filter[idxChannel] = new Filtfilt(13, a, b);
	}


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
	strcpy(filename_base, modelName_base);
	strcpy(filename_fn, modelName_fn);
	model_base = svm_load_model(filename_base);
	model_fn = svm_load_model(filename_fn);

	tmpfilename_base = strtok(filename_base, ".");
	strcat(tmpfilename_base, ".txt");
	tmpfilename_fn = strtok(filename_fn, ".");
	strcat(tmpfilename_fn, ".txt");

	fp_base = fopen(tmpfilename_base, "rt");
	fp_fn = fopen(tmpfilename_fn, "rt");

	for(int feat_index=0; feat_index < N_SelFeat; ++feat_index)
	{
		fscanf(fp_base, "%d", &featSel_base[feat_index]);
		fscanf(fp_fn, "%d", &featSel_fn[feat_index]);
	}

	fclose(fp_base);
	fclose(fp_fn);
	
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
		delete filter[i];
	}
	delete fftSpect;
	delete filter;
}

bool Classification::isFinished()
{
	return flagProcess;
}

void Classification::PutData(double *inSegment, int nCh, int nCol)
{
	for(int channel_index=0; channel_index<nCh; ++channel_index)
	{
		double *filtSegment;
		filtSegment = new double [nCol];

		//// filter data
		//filter[channel_index]->calculate(nCol, &inSegment[channel_index*nCol], filtSegment);

		// w/o filtering
		memcpy(filtSegment, &inSegment[channel_index*nCol], nCol*sizeof(double));
		
		
		// put data into queue
		if(headQue + nCol < QueueSize)
			memcpy(&queData[channel_index][headQue], filtSegment, nCol*sizeof(double));
		else
		{
			int diff = QueueSize-headQue;
			memcpy(&queData[channel_index][headQue], filtSegment, diff*sizeof(double));
			memcpy(&queData[channel_index][0], &filtSegment[diff], (nCol-diff)*sizeof(double));
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
	double tmpMax_base, tmpMin_base, tmpMax_fn, tmpMin_fn;

	diffHT = (headQue + QueueSize - tailQue_I) % QueueSize;
	while(diffHT >= DataProcSize)
	{
		// calculate FFT for DataProcSize (ex. 500ms)
		for(int channel_index=0; channel_index<nChannel; ++channel_index)
		{
			double inBuff[DataProcSize];

			// move moving window(data queue) to input buffer
			if(tailQue_I + DataProcSize <= QueueSize)
			{
				memcpy(inBuff, &queData[channel_index][tailQue_I], sizeof(double) * DataProcSize);
			}
			else
			{
				int nTail = QueueSize - tailQue_I;
				memcpy(inBuff, &queData[channel_index][tailQue_I], sizeof(double) * nTail);
				memcpy(&inBuff[nTail], &queData[channel_index][0], sizeof(double) * (DataProcSize-nTail));
			}

			// excute FFT
			memcpy(fftIn, inBuff, sizeof(double)* DataProcSize);
			fftw_execute(planFFT);

			// rearrange result FFT & append to spectrogram
			fftSpect[channel_index][0][idxSpect] = abs(fftOut[0]);
			for(int frequency_index=1;frequency_index<FFTSize/2;++frequency_index)
			{
				double Re = fftOut[frequency_index];
				double Im = fftOut[FFTSize-frequency_index];
				fftSpect[channel_index][frequency_index][idxSpect] = sqrt(Re*Re + Im*Im);
			}

		}

		++idxSpect;

		// when the data ready to make spectrogram (ex. 2s)
		if(idxSpect == SpectSize)
		{
			int cnt;
			int classResult;
			
		
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

							// reduce frequency dimension (divide into 5bands)
							for(int freqBin_index=freqRange[frequency_index][0];freqBin_index<=freqRange[frequency_index][1];++freqBin_index)
							{
								sum += fftSpect[channel_index][freqBin_index][time_index];
							}
							rdcSpect[channel_index][frequency_index][time_index] = sum / (freqRange[frequency_index][1] - freqRange[frequency_index][0] + 1);


							// moving average routine
							// accumulate spectral power along the time
							if(time_index == 0)
								accumSpect[channel_index][frequency_index][time_index] = rdcSpect[channel_index][frequency_index][time_index];
							else
								accumSpect[channel_index][frequency_index][time_index] = accumSpect[channel_index][frequency_index][time_index-1] + rdcSpect[channel_index][frequency_index][time_index];

							// 500ms moving average --> baseline
							if(time_index >= TS_500ms)
								accumSpect[channel_index][frequency_index][time_index] -= rdcSpect[channel_index][frequency_index][time_index-TS_500ms];
						}
					}
					for(int frequency_index=0;frequency_index<NoFreq;++frequency_index)
					{

						// making baseline and feature matrix

						//for(int time_index=0;time_index<NoTime;++time_index)
						//{
						//	double sum = 0;
						//	int seg_time = TS_500ms+(time_index*TS_100ms);
						//	for(int seg_index=seg_time; seg_index < seg_time+TS_100ms; ++seg_index)
						//	{
						//		 sum += rdcSpect[channel_index][frequency_index][seg_index];
						//	}
						//	featureMat[channel_index][frequency_index][time_index] = sum / TS_100ms;
						//}

						for(int time_index=0;time_index<NoTime;++time_index)
						{
							double sum = 0;
							int seg_time = TS_500ms+(6*TS_100ms);
							for(int seg_index=seg_time; seg_index < seg_time+(6*TS_100ms); ++seg_index)
							{
								sum += rdcSpect[channel_index][frequency_index][seg_index];
							}
							featureMat[channel_index][frequency_index][time_index] = sum / (6*TS_100ms);
						}

						baseMat[channel_index][frequency_index] = accumSpect[channel_index][frequency_index][TS_500ms-1] / TS_500ms;

						// normalize feature matrix by log10
						for(int time_index=0;time_index<NoTime;++time_index)
						{
							normMat[channel_index][frequency_index][time_index] = log10(featureMat[channel_index][frequency_index][time_index] / baseMat[channel_index][frequency_index]);
						}
					}
				}
				isFirst = false;
			}
			else	// second and after... -> calculate just for the shift
			{
				for(int channel_index=0;channel_index<nChannel;++channel_index)
				{
					for(int time_index=SpectSize-TS_100ms;time_index<SpectSize;++time_index)
					{
						for(int frequency_index=0;frequency_index<NoFreq;++frequency_index)
						{
							double sum = 0;

							// reduce frequency dimension
							for(int freqBin_index=freqRange[frequency_index][0];freqBin_index<=freqRange[frequency_index][1];++freqBin_index)
							{
								sum += fftSpect[channel_index][freqBin_index][time_index];
							}
							rdcSpect[channel_index][frequency_index][time_index] = sum / (freqRange[frequency_index][1] - freqRange[frequency_index][0] + 1);

							// accumulate spectral power along the time
							accumSpect[channel_index][frequency_index][time_index] = accumSpect[channel_index][frequency_index][time_index-1] + rdcSpect[channel_index][frequency_index][time_index];
							// 1000ms moving average 
							accumSpect[channel_index][frequency_index][time_index] -= rdcSpect[channel_index][frequency_index][time_index-TS_500ms];
						}
					}
					for(int frequency_index=0;frequency_index<NoFreq;++frequency_index)
					{
						// add feature matrix segment

						//double sum = 0;
						//int seg_time = TS_500ms+((NoTime-1)*TS_100ms);
						//for(int seg_index = seg_time; seg_index < seg_time+TS_100ms; ++seg_index)
						//{
						//	sum += rdcSpect[channel_index][frequency_index][seg_index];
						//}
						//featureMat[channel_index][frequency_index][NoTime-1] = sum / (TS_100ms);

						double sum = 0;
						int seg_time = TS_500ms+(6*TS_100ms);
						for(int seg_index = seg_time; seg_index < seg_time+(6*TS_100ms); ++seg_index)
						{
							sum += rdcSpect[channel_index][frequency_index][seg_index];
						}
						featureMat[channel_index][frequency_index][NoTime-1] = sum / (6*TS_100ms);

						// baseline matrix
						baseMat[channel_index][frequency_index] = accumSpect[channel_index][frequency_index][TS_500ms-1] / TS_500ms;
						// normalize the feature matrix segment
						normMat[channel_index][frequency_index][NoTime-1] = log10(featureMat[channel_index][frequency_index][NoTime-1] / baseMat[channel_index][frequency_index]);
					}
				}				
			}


			// linearize features (3d->1d, Ch x Time x Freq)
			nFeatures_tmp = nChannel * NoTime * NoFreq;
			nFeatures_base = N_SelFeat;
			nFeatures_fn = N_SelFeat;
			featureMat_tmp = new double [nFeatures_tmp];
			featureMat_base = new struct svm_node [nFeatures_base + 1];
			featureMat_fn = new struct svm_node [nFeatures_fn + 1];

			cnt = 0;
			for(int channel_index=0; channel_index<nChannel; ++channel_index)
			{
				for(int time_index=0; time_index<NoTime; ++time_index)
				{
					for(int freq_index=0; freq_index<NoFreq; ++freq_index)
					{
						featureMat_tmp[cnt++] = normMat[channel_index][freq_index][time_index];
					}
				}
			}


			// extract selected (high)features
			tmpMax_base = normMat[0][0][0];
			tmpMin_base = tmpMax_base;
			tmpMax_fn = tmpMax_base;
			tmpMin_fn = tmpMax_base;
			for(int feature_index=0; feature_index < N_SelFeat; ++feature_index)
			{
				featureMat_base[feature_index].index = feature_index+1;
				featureMat_base[feature_index].value = featureMat_tmp[featSel_base[feature_index]];
				featureMat_fn[feature_index].index = feature_index+1;
				featureMat_fn[feature_index].value = featureMat_tmp[featSel_fn[feature_index]];
				if(tmpMax_base < featureMat_base[feature_index].value) tmpMax_base = featureMat_base[feature_index].value;
				if(tmpMin_base > featureMat_base[feature_index].value) tmpMin_base = featureMat_base[feature_index].value;
				if(tmpMax_fn < featureMat_fn[feature_index].value) tmpMax_fn = featureMat_fn[feature_index].value;
				if(tmpMin_fn > featureMat_fn[feature_index].value) tmpMin_fn = featureMat_fn[feature_index].value;
			}
			featureMat_base[nFeatures_base].value = -1;
			featureMat_fn[nFeatures_base].value = -1;
			tmpMax_base -= tmpMin_base;
			tmpMax_fn -= tmpMin_fn;


			// feature normalization (to -1~1) : ((feat - min) / max) * 2 - 1 
			for(int feature_index=0; feature_index < N_SelFeat; ++feature_index)
			{
				featureMat_base[feature_index].value = (featureMat_base[feature_index].value - tmpMin_base) / tmpMax_base * 2 - 1;
				featureMat_fn[feature_index].value = (featureMat_fn[feature_index].value - tmpMin_fn) / tmpMax_fn * 2 - 1;
			}

			// classification
			// test if the block is baseline or not
			classResult = (int)svm_predict(model_base, featureMat_base);
			//classResult = 1;
			if(classResult)
			{
				// differentiate face/number
				classResult = (int)svm_predict(model_fn, featureMat_fn);
			}


			//// result voting
			//++labelCnt[classResult];
			//if(cntLabels < NoLabel-1)
			//{// initiallize vote
			//	recentLabels[cntLabels++] = classResult;
			//	iResult = 0;
			//}
			//else
			//{
			//	int labelMax = 0;

			//	recentLabels[cntLabels] = classResult;
			//	
			//	// voting system (more than half, default : neutral) 
			//	if(labelCnt[labelMax] < labelCnt[1])
			//		labelMax = 1;
			//	if(labelCnt[labelMax] < labelCnt[2])
			//		labelMax = 2;
			//	if(labelCnt[labelMax] < NoLabel/2)
			//		labelMax = 0;
			//	iResult = labelMax;

			//	// shift labels for next use
			//	labelCnt[recentLabels[0]]--;
			//	for(int label_index=0;label_index<NoLabel-1;++label_index)
			//	{
			//		recentLabels[label_index] = recentLabels[label_index+1];
			//	}
			//}

			iResult = classResult;

			delete featureMat_base;
			
			
			// shift 100ms
			for(int channel_index=0;channel_index<nChannel;++channel_index)
			{
				for(int frequency_index=0;frequency_index<(int)(FFTSize/2);++frequency_index)
				{
					for(int time_index=0;time_index<SpectSize-TS_100ms;++time_index)
					{
						fftSpect[channel_index][frequency_index][time_index] = fftSpect[channel_index][frequency_index][time_index+TS_100ms];
					}
				}
				for(int frequency_index=0;frequency_index<NoFreq;++frequency_index)
				{
					for(int time_index=0;time_index<SpectSize-TS_100ms;++time_index)
					{
						rdcSpect[channel_index][frequency_index][time_index] = rdcSpect[channel_index][frequency_index][time_index+TS_100ms];
						accumSpect[channel_index][frequency_index][time_index] = accumSpect[channel_index][frequency_index][time_index+TS_100ms];
					}
					for(int time_index=0;time_index<NoTime-1;++time_index)
					{
						featureMat[channel_index][frequency_index][time_index] = featureMat[channel_index][frequency_index][time_index+1];
						normMat[channel_index][frequency_index][time_index] = normMat[channel_index][frequency_index][time_index+1];
					}
				}
			}
			idxSpect -= TS_100ms;
		}

		// shift moving window
		tailQue += T_50ms;
		if(tailQue_I >= QueueSize)
			while((tailQue -= QueueSize) >= QueueSize);
		diffHT = (headQue + QueueSize - tailQue_I) % QueueSize;
	}

	flagProcess = 1;

	return iResult;
}



