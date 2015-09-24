#pragma once
#include<Windows.h>
#include "fftw3.h"
#include "svm.h"

#define SamplingRate	256
#define DataProcSize	((int)(SamplingRate * 0.5))			// 128 |256sps (500ms)
//#define T_20ms			(SamplingRate * 0.02)				// 5.12 |256sps (50ms)
//#define TS_1000ms		((int)(SamplingRate / T_20ms))		// 50 |256sps (1s)
//#define TS_100ms		((int)(TS_1000ms * 0.1))			// 5 |256sps (150ms)
//#define TS_500ms		((int)(TS_1000ms * 0.5))			// 25 |256sps (500ms)
#define T_50ms			(SamplingRate * 0.05)				// 12.5 |250sps (50ms)
#define TS_1000ms		((int)(SamplingRate / T_50ms))		// 20 |250sps (1s)
#define TS_100ms		((int)(TS_1000ms * 0.1))			// 2 |250sps (150ms)
#define TS_500ms		((int)(TS_1000ms * 0.5))			// 10 |250sps (500ms)
#define QueueSize		1000
#define NoChannel		16
#define NoFreq			5
#define NoTime			9
#define NoLabel			7
#define FFTSize			256
#define SpectSize		40
#define N_SelFeat		100
#define tailQue_I		((int)tailQue)

typedef fftw_plan (*pPlan1d)(int, double*, double*, fftw_r2r_kind, unsigned);
typedef void (*pExecute)(const fftw_plan);

const int freqRange[5][2] = {{0,2},{3,6},{7,10},{11,28},{29,49}};

class Classification
{
public:
	Classification(void);
	// initiallization with SVM coefficients
	Classification(char*, char*);
	~Classification(void);

	// process functions
	void PutData(double *, int, int);
	int ProcessData();
	bool isFinished();


private:
	bool flagProcess;
	int iResult;

	// SVM
	struct svm_model *model_base, *model_fn;
	struct svm_node *featureMat_base, *featureMat_fn;
	double *featureMat_tmp;
	int featSel_base[N_SelFeat], featSel_fn[N_SelFeat];
	int nFeatures_tmp, nFeatures_base, nFeatures_fn;
//	struct svm_coeff *idx_base, *idx_fn;

	// data queue
	double queData[NoChannel][QueueSize];
	int nChannel;
	int headQue;
	double tailQue;

	// variables for DLL load
	HINSTANCE hFFTW;
	pPlan1d fftw_plan_r2r_1d;
	pExecute fftw_execute;
	fftw_plan planFFT;

	// variables for FFT
	double fftOut[FFTSize];
	double fftIn[FFTSize];

	// Spectrogram
	double ***fftSpect;
	int idxSpect;
	double baseMat[NoChannel][NoFreq];
	double featureMat[NoChannel][NoFreq][NoTime];
	double normMat[NoChannel][NoFreq][NoTime];
	double accumSpect[NoChannel][NoFreq][SpectSize];
	double rdcSpect[NoChannel][NoFreq][SpectSize];
	bool isFirst;

	// vote algorithm
	int recentLabels[NoLabel];
	int cntLabels;
	int labelCnt[3];
};

