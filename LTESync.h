#ifndef _LTESYNC_H_
#define _LTESYNC_H_

#include <fstream>
using namespace std;

#define pi 3.14159265

extern int PSS_Peak_idx;
extern int sf;

void ifft(float *inputI, float *inputQ, int length_fft, float *outputI, float *outputQ);
void fft(float *inputI, float *inputQ, int scno, float *outputI, float *outputQ);

void LTESync_init();
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

int mod(float a, int b);
int mod(int a, int b);

#endif