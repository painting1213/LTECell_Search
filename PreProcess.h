#ifndef _PREPROCESS_H_
#define _PREPROCESS_H_

#include <fstream>
using namespace std;

// 打开文件，文件有两列数据：I 路， Q 路
//	sscache = 对数据1/16重采样，使用半带滤波器HBF或CIC---------滤波, 30.72e6*0.001=30.72e3, 30.72e3/16=1920
//	sscache_f = 对sscache进行滤波(卷积)，给定filter_coffe系数，得到的数据与sscache长度相同-------卷积，滤波
//  sscache_d = 对卷积后的结果进行1抽取，使用半带滤波器HBF-------滤波
//struct
struct Signal_Short
{
    short I_Data; //2 bytes
    short Q_Data;
};

struct Signal_Float
{
    float I_Data; //4 bytes
    float Q_Data;
};

const int NumSamples_total = 184320;      //6ms， 30.72e6 * 0.06 = 184320
const int NumSamples_subframe = 30720;    //read Numamples numbers every 1ms

extern float I_total[NumSamples_total];               //get the total data
extern float Q_total[NumSamples_total];               //get the total data

extern float ss_cs_I[NumSamples_subframe];            //get 30720 numbers
extern float ss_cs_Q[NumSamples_subframe];            //get 30720 numbers

extern float ss_cs_I_16Dec[NumSamples_subframe / 16];   //after the HBF*4
extern float ss_cs_Q_16Dec[NumSamples_subframe / 16];   //after the HBF*4

extern float ss_cs_I_32Dec[NumSamples_subframe / 32]; // the data for following process
extern float ss_cs_Q_32Dec[NumSamples_subframe / 32]; // the data for following process

void PreProcess_init();
void WriteData();

void Get_Short_Data(ifstream &infile);
void Get_Float_Data(ifstream &infile);
//void GetData(ifstream &infile);                  //&infile

void Resample(float *ss_cs_I, float *ss_cs_Q, float *ss_cs_I_16Dec, float *ss_cs_Q_16Dec);
void Halfband_filter(float *input, int inlength, const float *filtercoef, int filength, float *output);
void Decimation(float *input, int inlength, int Deci, float *output);
void Conv(float *x, int x_length, const float *h, int h_length,float *output);
void Conv_same(float *x, int x_length, const float *h, const int h_length, float *output);

#endif