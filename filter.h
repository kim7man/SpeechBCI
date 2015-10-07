#pragma once
#include<Windows.h>


class Filtfilt
{
public:
	Filtfilt(void);
	// initiallization with filter coefficients
	Filtfilt(int n, const double*, const double*);
	~Filtfilt(void);

	void calculate(int, double *, double *);
	double calculate(double);


private:
	double *numer;
	double *denom;

	int nCoeff;

	double *x_delay;
	double *y_out;

	double MultiCoeff(double);

};
