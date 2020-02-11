#ifndef _PREPROCESS_H_
#define _PREPROCESS_H_

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include "Head.h"
using namespace std;

// 打开文件，文件有两列数据：I 路， Q 路
//	sscache = 对数据1/16重采样，使用半带滤波器HBF或CIC---------滤波, 30.72e6*0.001=30.72e3, 30.72e3/16=1920
//	sscache_f = 对sscache进行滤波(卷积)，给定filter_coffe系数，得到的数据与sscache长度相同-------卷积，滤波
//  sscache_d = 对卷积后的结果进行1抽取，使用半带滤波器HBF-------滤波


const int NumSamples_total = 184320;      //6ms， 30.72e6 * 0.06 = 184320
const int NumSamples_subframe = 30720;    //read Numamples numbers every 1ms

const int hbfiltercoef2length = 7;      //const
const int hbfiltercoef1length = 11;
const int decfilterlength = 17;
const int filter_ds_length = 11;

const float hbfiltercoef2[7] = {1.0*-1058/16384, 0, 1.0*9250/16384, 1.0, 1.0*9250/16384, 0, 1.0*-1058/16384};
const float hbfiltercoef1[11] = {1.0*242/16384, 0, 1.0*-1736/16384, 0, 1.0*9686/16384, 1.0, 1.0*9686/16384, 0,
                                1.0*-1736/16384, 0, 1.0*242/18384};
const float decfiltercoef[17] = {-80.0/15394,-146.0/15394,287.0/15394,742.0/15394,-589.0/15394,-2449.0/15394, 
                                877.0/15394,10045.0/15394,1.0,10045.0/15394,877.0/15394,-2449.0/15394,-589.0/15394, 
                                742.0/15394,287.0/15394,-146.0/15394,-80.0/15394};
const float filter_ds[]={0.125, 0, -0.25, 0, 0.625, 1, 0.625, 0, -0.25, 0, 0.125};

float I_total[NumSamples_total];               //get the total data
float Q_total[NumSamples_total];               //get the total data

float ss_cs_I[NumSamples_subframe];            //get 30720 numbers
float ss_cs_Q[NumSamples_subframe];            //get 30720 numbers

float ss_cs_I_Dec[NumSamples_subframe / 16];   //after the HBF*4
float ss_cs_Q_Dec[NumSamples_subframe / 16];   //after the HBF*4

float ss_cs_I_process[NumSamples_subframe / 32]; // the data for following process
float ss_cs_Q_process[NumSamples_subframe / 32]; // the data for following process

void GetData(ifstream &infile);                  //&infile
void PreProcess_init(ifstream &infile);

void Resample(float *ss_cs_I, float *ss_cs_Q, float *ss_cs_I_dec, float *ss_cs_Q_dec);
void halfband_filter(float *input, int inlength, const float *filtercoef, int filength, float *output);
void Decimation(float *input, int inlength, int Deci, float *output);
void conv(float *x, int x_length, const float *h, int h_length,float *output);
void conv_same(float *x, int x_length, const float *h, const int h_length, float *output);

#endif