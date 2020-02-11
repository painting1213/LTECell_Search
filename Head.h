#ifndef _HEAD_H_
#define _HEAD_H_

#include <iostream>
#include<cmath>
using namespace std;

#define pi 3.14159265

int mod(float a, int b);
int mod(int a, int b);
void return_zero(float *x, int len);
void FIR_filter(double *input, int inputlength, short *filtercoef, int filterlength, double *output);

#endif