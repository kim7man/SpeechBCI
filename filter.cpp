#include<memory.h>

#include "filter.h"

Filtfilt::Filtfilt(void)
{
	numer = NULL;
	denom = NULL;
}

Filtfilt::~Filtfilt(void)
{
	delete numer;
	delete denom;
	delete x_delay;
}

Filtfilt::Filtfilt(int n, const double *a, const double *b)
{
	nCoeff = n;

	numer = new double [nCoeff];
	denom = new double [nCoeff];
	x_delay = new double[nCoeff];
	y_out = new double[nCoeff];

	memcpy(numer, b, sizeof(double) * nCoeff);
	memcpy(denom, a, sizeof(double) * nCoeff);
	memset(x_delay, 0, sizeof(double) * nCoeff);
	memset(y_out, 0, sizeof(double) * nCoeff);

}


void Filtfilt::calculate(int n, double *src, double *ret)
{
	int nData = n;
	double *data, *result;

	data = new double [nData];
	result = new double [nData];
	memcpy(data, src, sizeof(double) * nData);

	for(int idxData=0; idxData<nData; ++idxData)
	{
		result[idxData] = MultiCoeff(data[idxData]);
	}
	
	memcpy(ret, result, sizeof(double) * nData);

	delete data;
	delete result;
}

double Filtfilt::MultiCoeff(double newData)
{

	// shift the old samples
	for(int idxCoeff=nCoeff-1; idxCoeff>0; --idxCoeff)
	{
		x_delay[idxCoeff] = x_delay[idxCoeff-1];
		y_out[idxCoeff] = y_out[idxCoeff-1];
	}

	// put new data
	x_delay[0] = newData;

	// multiply coefficients
	y_out[0] = 0;
	for(int idxCoeff=0; idxCoeff<nCoeff; ++idxCoeff)
	{
		y_out[0] += x_delay[idxCoeff] * numer[idxCoeff];
	}
	for(int idxCoeff=1; idxCoeff<nCoeff; ++idxCoeff)
	{
		y_out[0] -= y_out[idxCoeff] * denom[idxCoeff];
	}
	y_out[0] /= denom[0];
	
	return y_out[0];
}
