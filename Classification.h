#pragma once
#include<Windows.h>
#include "fftw3.h"
#include "svm.h"

#define SamplingRate	250
#define DataProcSize	((int)(SamplingRate / 2))	// 125 |250sps
#define T_50ms			((int)(SamplingRate * 0.05))	// 12 |250sps
#define TS_1000ms		((int)(SamplingRate / T_50ms))	// 20 |250sps
#define TS_150ms		((int)(SamplingRate * 0.15 / T_50ms))	// 3 |250sps
#define QueueSize		300
#define NoChannel		30
#define NoFreq			5
#define NoTime			7
#define NoLabel			7
#define FFTSize			256
#define SpectSize		40

typedef fftw_plan (*pPlan1d)(int, double*, double*, fftw_r2r_kind, unsigned);
typedef void (*pExecute)(const fftw_plan);

const int freqRange[5][2] = {{0,3},{4,7},{8,11},{12,30},{30,50}};

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
	int nFeatures_base, nFeatures_fn;
	struct svm_coeff *idx_base, *idx_fn;

	// data queue
	double queData[NoChannel][QueueSize];
	int nChannel;
	int headQue;
	int tailQue;

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

