#include "Head.h"

//make all x[len] = 0
void return_zero(float *x, int len)
{
    for (int i = 0; i < len; i++)
    {
        x[i] = 0;
    }
}

//mod(float, int)
int mod(float a, int b)
{
	int temp = int(a);
	return temp % b;
}

//mod(int, int)
int mod(int a, int b)
{
	return a % b;
}

// FIR Filter
void FIR_filter(double *input, int inputlength, short *filtercoef, int filterlength, double *output)
{
	int i = 0,j = 0,k=0;
	double temp; 
    double state[filterlength+1];

	for (k = 0; k < inputlength; k++) 
	{ 
		state[0] = input[k];
		for (i = 0, temp = 0; i < filterlength; i++)
			temp += filtercoef[i] * state[i];
		    output[k] = temp;
		for (j = filterlength - 1; j > -1 ; j--)
		state[j+1] = state[j];
	}
}