#ifndef _LTESYNC_H_
#define _LTESYNC_H_

#include <iostream>
#include <fstream>
using namespace std;

#define pi 3.14159265

const int scno = 64;
const int Length_FFT = 64;
const int Length_PSS = 62;
const int Nums_subframes = 6;
const int Nums_per_subframe = 960; //每个子帧1ms,每一个子帧的样本点数目=30.72M*1/1000*32 = 960

//const short subframe_gap = 6;
//const short subframe_cache = 6;

int PSS_Peak_idx = 0;
int sf = 0;

//PSS frequency domain data
float PSSfreq_I[62];     //after genPSS(Nid_2)
float PSSfreq_Q[62];

float SSS_du_0[62];      //after genSSS(Nid_1, Nid_2)
float SSS_du_5[62];

//float corr[5760][3] = {0};

float pss_preprocess_I[5760];
float pss_preprocess_Q[5760]; 

float sss_preprocess_I[5760 * 2];
float sss_preprocess_Q[5760 * 2];


float corr[5696][3];

void ifft(float *inputI, float *inputQ, int length_fft, float *outputI, float *outputQ);
void fft(float *inputI, float *inputQ, int scno, float *outputI, float *outputQ);

void LTESync_init(ifstream &infile_pss, ifstream &infile_sss);
void get_preprocss_data(ifstream &infile, float *preprocess_I, float *preprocess_Q);                            //ss_cs_I[5760], ss_cs_Q[5760]
void swapPSS(float *PSSfreq_I, float *PSSfreq_Q, float *PSSfreq_tmpI, float *PSSfreq_tmpQ);
//void GetPSSfreq(float (*PSSfreq_I)[62], float (*PSSfreq_Q)[62], int N_id_2,float *pssfreq_I, float *pssfreq_Q);

float calc_Sequence_Modulus(float *ss_cs_I, float *ss_cs_Q, int samp_idx);
float calc_complex_Multi(float *PSStime_I, float *PSStime_Q, float *ss_cs_I, float *ss_cs_Q, int samp_idx);

//generate PSS, SSS
void genPSS(int Nid_2);                   //get PSSfreq_I[62], PSSfreq_Q[62] 
void genSSS(int Nid_1, int Nid_2);        //get SSS_du_0[62],  SSS_du_5[62]

//find SSS_signal, ss_cs[5760]
int findPSS(float *ss_cs_I, float *ss_cs_Q);                               //return Nid_2
int findSSS(float *ss_cs_I, float *ss_cs_Q, int Nid_2, int pss_peak_idx);  //return Nid_1

int ident_Nid_2(float (*corr)[3]);
int Set_CellID(int Nid_1, int Nid_2);

#endif