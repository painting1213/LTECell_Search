#include <iostream>
#include <string.h>
#include <cmath>

#include "LTESync.h"

//Define gloabal variables
int PSS_Peak_idx = 0;
int sf;

//const variables
const int scno = 64;
const int Length_FFT = 64;
const int Length_PSS = 62;
const int Nums_subframes = 6;
const int Nums_per_subframe = 960; //每个子帧1ms,每一个子帧的样本点数目=30.72M*1/1000*32 = 960
//const short subframe_gap = 6;
//const short subframe_cache = 6;

//PSS frequency domain data
float PSSfreq_I[62];     //after genPSS(Nid_2)
float PSSfreq_Q[62];
float corr[5696][3];

//SSS data
float SSS_du_0[62];      //after genSSS(Nid_1, Nid_2)
float SSS_du_5[62];

//process data
float pss_preprocess_I[5760];
float pss_preprocess_Q[5760]; 

float sss_preprocess_I[5760 * 2];
float sss_preprocess_Q[5760 * 2];

//LTE Sync Init
void LTESync_init()
{
	ifstream infile_pss("usrp_process_pss.dat", ios::in | ios::binary);
	ifstream infile_sss("usrp_process_sss.dat", ios::in | ios::binary);

	get_preprocss_data(infile_pss, pss_preprocess_I, pss_preprocess_Q); //get pss_preprocess[5760]
	get_preprocss_data(infile_sss, sss_preprocess_I, sss_preprocess_Q); //get sss_preprocess[5760 * 2]
}

//From File to get preprocess_I, preprocess_I
void get_preprocss_data(ifstream &infile, float *preprocess_I, float *preprocess_Q)
{
	struct Signal_float
	{
		float I_data; // 4 bytes
		float Q_data; // 4 bytes
	};

	Signal_float signal;
	int readBytes = 0;

	if (!infile)
	{
		std::cout << "Error, can not open the file \n"
				  << endl;
		return;
	}

	while (infile.read((char *)&signal, sizeof(signal)))
	{
		if (readBytes % 8 == 0 && readBytes != 0)
		{
			preprocess_I[readBytes / 8] = signal.I_data;
			preprocess_Q[readBytes / 8] = signal.Q_data;
		}
		readBytes += infile.gcount(); //count the read bytes
	}
	//cout << "the numbers of total data is: " << readBytes / 8 << endl;
	infile.close();
}

//output = ifft(input), the length of sequence is scno
void ifft(float *inputI, float *inputQ, int scno, float *ouputI, float *ouputQ)
{
	int n = 0;
	int k = 0;

	//float i = 0;
	//float q = 0;
	float sumi = 0;
	float sumq = 0;
	float expi = 0;
	float expq = 0;

	for (n = 0; n < scno; n++)
	{
		sumi = 0;
		sumq = 0;

		for (k = 0; k < scno; k++)
		{
			expi = cos(1 * 2 * pi * n * k / scno);
			expq = sin(1 * 2 * pi * n * k / scno);
			sumi = sumi + inputI[k] * expi - inputQ[k] * expq;
			sumq = sumq + inputI[k] * expq + inputQ[k] * expi;
		}

		ouputI[n] = sumi / (1.0 * scno);
		ouputQ[n] = sumq / (-1.0 * scno);
	}

	return;
}

/***********Perform FFT operation********/
void fft(float *inputI, float *inputQ, int scno, float *outputI, float *outputQ)
{
	int n = 0;
	int k = 0;

	//double i = 0;
	//double q = 0;
	double sumi = 0;
	double sumq = 0;
	double expi = 0;
	double expq = 0;

	for (k = 0; k < scno; k++)
	{
		sumi = 0;
		sumq = 0;

		for (n = 0; n < scno; n++)
		{
			expi = cos(1 * 2 * pi * k * n / scno);
			expq = sin(-1 * 2 * pi * k * n / scno);

			sumi = sumi + inputI[n] * expi - inputQ[n] * expq;
			sumq = sumq + inputI[n] * expq + inputQ[n] * expi;
		}
		outputI[k] = sumi;
		outputQ[k] = sumq;
	}
	return;
}

/*******swap PSSfreq_I and PSSfreq_Q ,the results are PSSfreq_tmpI and PSSfreq_tmpQ*******/
// PSSfreq_I[62], PSSfreq_Q[62],  and PSSfreq_tmpI[64], PSSfreq_tmpQ[64]
void swapPSS(float *PSSfreq_I, float *PSSfreq_Q, float *PSSfreq_tmpI, float *PSSfreq_tmpQ)
{
	int i = 0;
	// int length = 64;

	// float *psstmpI = new float[length]; //delete
	// float *psstmpQ = new float[length]; //delete

	// memset(psstmpI, 0, length); // make all tmpI = 0
	// memset(psstmpQ, 0, length); // make all tmpQ = 0

	for (i = 1; i < 32; i++)
	{
		PSSfreq_tmpI[i] = PSSfreq_I[30 + i]; //tmp(1:31) = PSSfreq(31:61), begin at tmp[1]
		PSSfreq_tmpQ[i] = PSSfreq_Q[30 + i];

		PSSfreq_tmpI[i + 32] = PSSfreq_I[i - 1]; //tmp(33:63) = PSSfreq(0:30), begin at tmp[33]
		PSSfreq_tmpQ[i + 32] = PSSfreq_Q[i - 1];
	}

	// for (i = 0; i < length; i++)
	// {
	// 	PSSfreq_tmpI[i] = psstmpI[i];
	// 	PSSfreq_tmpQ[i] = psstmpQ[i];
	// }

	// delete[] psstmpI;
	// delete[] psstmpQ;
}

/**************calculate the modulus of complex sequence*************/
// the beginning is samp_idx, ss_cs[samp_idx]
float calc_Sequence_Modulus(float *ss_cs_I, float *ss_cs_Q, int samp_idx)
{
	float Modulus = 0;
	for (int i = samp_idx; i < samp_idx + Length_FFT; i++) //the length of sequence is 64
	{
		Modulus = Modulus + ss_cs_I[i] * ss_cs_I[i] + ss_cs_Q[i] * ss_cs_Q[i];
	}

	return Modulus;
}

/**************calculate the multi of two complex sequences*************/
// the beginning is samp_idx, ss_cs[samp_idx]
float calc_complex_Multi(float *PSStime_I, float *PSStime_Q, float *ss_cs_I, float *ss_cs_Q, int samp_idx)
{
	float Multi = 0;
	float tmpI = 0;
	float tmpQ = 0;

	for (int i = samp_idx, j = 0; i < samp_idx + Length_FFT, j < 64; i++, j++) //the length of sequence is 64
	{
		//Multi = Multi + pow(fabs(PSStime_I[j] * ss_cs_I[i] + PSStime_Q[j] * ss_cs_Q[i]), 2);
		tmpI = tmpI + PSStime_I[j] * ss_cs_I[i] - PSStime_Q[j] * ss_cs_Q[i];
		tmpQ = tmpQ + PSStime_I[j] * ss_cs_Q[i] + PSStime_Q[j] * ss_cs_I[i];
	}

	Multi = (tmpI * tmpI + tmpQ * tmpQ) * 1;
	//整数绝对值    int c = abs(b-c);
	//浮点数绝对值  double d = fabs(b-c)
	return Multi;
}

//get freq[62] from PSSfreq[3][62]
/*void GetPSSfreq(float (*PSSfreq_I)[62], float (*PSSfreq_Q)[62], int N_id_2, float *freq_I, float *freq_Q)
{
	if (N_id_2 >= 0 && N_id_2 <= 2)
	{
		for (int i = 0; i < 62; i++)
		{
			freq_I[i] = PSSfreq_I[N_id_2][i];
			freq_Q[i] = PSSfreq_Q[N_id_2][i];
		}
	}
	else
	{
		cout << "N_id_2 is wrong, please check!" << endl;
	}
}*/

/********************Generate PSSfreq_I[62], PSSfreq_Q[62] Sequence***********************/
//genPSS(Nid_2)
void genPSS(int Nid_2)
{
	int root_idx[] = {25, 29, 34}; //root_idx
	int u = root_idx[Nid_2];

	//validate Nid_2
	if (!(Nid_2 >= 0 && Nid_2 <= 2))
	{
		std::cout << "Error: Invalid Nid_2: " << Nid_2 << endl;
		return;
	}

	for (int i = 0; i < 31; i++)
	{
		PSSfreq_I[i] = cos(-1 * pi * u * i * (i + 1) / 63.0);
		PSSfreq_Q[i] = sin(1 * pi * u * i * (i + 1) / 63.0);
	}
	for (int i = 31; i < 62; i++)
	{
		PSSfreq_I[i] = cos(-1 * pi * u * (i + 2) * (i + 1) / 63.0);
		PSSfreq_Q[i] = sin(1 * pi * u * (i + 2) * (i + 1) / 63.0);
	}
}

/********************Generate the SSS_du_0[62], SSS_du_5[62]**********************/
// genSSS(Nid_1, Nid_2)
void genSSS(int Nid_1, int Nid_2)
{
	float q = 0;
	float m0 = 0;
	float m1 = 0;
	float q_prime = 0;
	float m_prime = 0;

	int s0_m0_idx = 0; //s0_m0[idx]
	int s1_m1_idx = 0; //s1_m1[idx]

	int c0_idx = 0; //c0[idx]
	int c1_idx = 0; //c1[idx]

	int z1_m0_idx = 0; //z1_m1[idx]
	int z1_m1_idx = 0; //z1_m1[idx]

	int x_s_tilda[31] = {0, 0, 0, 0, 1}; // [0]=0, [1]=0, [2]=0, [3]=0, [4]=1
	int s_tilda[31];					 //not initialize

	int x_c_tilda[31] = {0, 0, 0, 0, 1}; // [0]=0, [1]=0, [2]=0, [3]=0, [4]=1
	int c_tilda[31];					 //not initialize

	int x_z_tilda[31] = {0, 0, 0, 0, 1}; // [0]=0, [1]=0, [2]=0, [3]=0, [4]=1
	int z_tilda[31];

	int s0_m0[31];
	int s1_m1[31];

	int c0[31];
	int c1[31];

	int z1_m0[31];
	int z1_m1[31];

	//validate Nid_1 and Nid_2
	if (!(Nid_1 >= 0 && Nid_1 <= 167))
	{
		std::cout << "Error, Invalid Nid_1: " << Nid_1 << endl;
		return;
	}

	if (!(Nid_2 >= 0 && Nid_2 <= 2))
	{
		std::cout << "Error, Invalid Nid_2: " << Nid_2 << endl;
		return;
	}

	// generate m0 and m1
	q_prime = floor((1.0 * Nid_1) / 30); //floor(x)<=x
	q = floor(1.0 * (Nid_1 + q_prime * (q_prime + 1) / 2) / 30);
	m_prime = Nid_1 + q * (q + 1) / 2;
	m0 = mod(m_prime, 31);
	m1 = mod((m0 + floor(m_prime / 31) + 1), 31);

	// generate s_tilda
	for (int i = 0; i < 26; i++)
	{
		x_s_tilda[i + 5] = mod((x_s_tilda[i + 2] + x_s_tilda[i]), 2);
	}
	for (int j = 0; j < 31; j++)
	{
		s_tilda[j] = 1 - 2 * x_s_tilda[j];
	}

	// generate c_tilda
	for (int i = 0; i < 26; i++)
	{
		x_c_tilda[i + 5] = mod((x_c_tilda[i + 3] + x_c_tilda[i]), 2);
	}
	for (int j = 0; j < 31; j++)
	{
		c_tilda[j] = 1 - 2 * x_c_tilda[j];
	}

	// generate z_tilda
	for (int i = 0; i < 26; i++)
	{
		x_z_tilda[i + 5] = mod((x_z_tilda[i + 4] + x_z_tilda[i + 2] + x_z_tilda[i + 1] + x_z_tilda[i]), 2);
	}
	for (int j = 0; j < 31; j++)
	{
		z_tilda[j] = 1 - 2 * x_z_tilda[j];
	}

	// generate s0_m0[31], s1_m1[31]
	for (int n = 0; n < 31; n++)
	{
		s0_m0_idx = mod((n + m0), 31);
		s1_m1_idx = mod((n + m1), 31);

		s0_m0[n] = s_tilda[s0_m0_idx];
		s1_m1[n] = s_tilda[s1_m1_idx];
	}

	// generate c0,c1
	for (int n = 0; n < 31; n++)
	{
		c0_idx = mod((n + Nid_2), 31);
		c1_idx = mod((n + Nid_2 + 3), 31);

		c0[n] = c_tilda[c0_idx];
		c1[n] = c_tilda[c1_idx];
	}

	// generate z1_m0, z1_m1
	for (int n = 0; n < 31; n++)
	{
		z1_m0_idx = mod(n + (mod(m0, 8)), 31);
		z1_m1_idx = mod(n + (mod(m1, 8)), 31);

		z1_m0[n] = z_tilda[z1_m0_idx];
		z1_m1[n] = z_tilda[z1_m1_idx];
	}

	// generate SSS
	// SSS_du_0[62], SSS_du_5[62]
	for (int n = 0; n < 31; n++)
	{
		SSS_du_0[2 * n] = s0_m0[n] * c0[n];
		SSS_du_5[2 * n] = s1_m1[n] * c0[n];

		SSS_du_0[2 * n + 1] = s1_m1[n] * c1[n] * z1_m0[n];
		SSS_du_5[2 * n + 1] = s0_m0[n] * c1[n] * z1_m1[n];
	}
}

//ss_cs_I[5760], ss_cs_Q[5760], return Nid_2
int findPSS(float *ss_cs_I, float *ss_cs_Q)
{
	int Nid_2;
	//int idx = Nums_subframes * Nums_per_subframe - Length_FFT -1;
	int idx = 5760 - 64 - 1;
	float xx = 0, xy = 0, yy = 0, nv = 0;

	//float corr[5696][3];
	//memset(corr, 0, 5760*3);

	//temp[64] of PSSfreq_I, PSSfreq_Q
	float *PSSfreq_tmpI = new float[64];
	float *PSSfreq_tmpQ = new float[64];

	memset(PSSfreq_tmpI, 0, 64); // make PSSfreq_tmpI[0-63] = 0
	memset(PSSfreq_tmpQ, 0, 64);

	//PSS time domain data
	float *PSStime_I = new float[64]; //64 points, 0 ~ 63
	float *PSStime_Q = new float[64];

	memset(PSStime_I, 0, 64); // make PSStime_I[0-63] = 0
	memset(PSStime_Q, 0, 64);

	for (Nid_2 = 0; Nid_2 < 3; Nid_2++)
	{
		genPSS(Nid_2);														// get PSSfreq_I[62], get PSSfreq_Q[62]
		swapPSS(PSSfreq_I, PSSfreq_Q, PSSfreq_tmpI, PSSfreq_tmpQ);			// get PSSfreq_temp[64]
		ifft(PSSfreq_tmpI, PSSfreq_tmpQ, Length_FFT, PSStime_I, PSStime_Q); // get PSS_time[64]

		for (int i = 0; i < Length_FFT; i++)
		{
			PSStime_Q[i] = -1.0 * PSStime_Q[i]; // conj(a + bi) = a - bi
		}

		xx = calc_Sequence_Modulus(PSStime_I, PSStime_Q, 0); //pss_td'*pss_td

		for (int samp_idx = 0; samp_idx < idx; samp_idx++)
		{
			yy = calc_Sequence_Modulus(pss_preprocess_I, pss_preprocess_Q, samp_idx);
			xy = calc_complex_Multi(PSStime_I, PSStime_Q, pss_preprocess_I, pss_preprocess_Q, samp_idx);
			nv = yy - xy / xx;

			corr[samp_idx][Nid_2] = xy / nv;
		}
	}

	Nid_2 = ident_Nid_2(corr);

	delete[] PSSfreq_tmpI;
	delete[] PSSfreq_tmpQ;

	delete[] PSStime_I;
	delete[] PSStime_Q;

	// ofstream outfile("data_corr_5696_3.txt", ios::out);
	// for (int ii = 0; ii < 5696; ii++)
	// {
	// 	outfile << corr[ii][0] << " " << corr[ii][1] << " " <<corr[ii][2] << endl;
	// }
	// outfile.close();

	return Nid_2;
}

int ident_Nid_2(float (*corr)[3])
{
	//  分别找corr三列数组中每一列最大的绝对值val0 val1 val2，和该最大值的索引值idx0 idx1 idx2,再求三列中的最大值
	//  Nid_2 = 最大值所在列值i
	//  下面求peak_idx, 如果idx(1) < 128*4,peak_idx=idx(i)+n_samps_per_subframe*5,否则等于idx(i)

	int i, j;
	int Nid_2;
	// int idx_max = 0;

	float val[3]; //每一列最大的绝对值val0 val1 val2
	int idx[3];   //该最大值的索引值idx0 idx1 idx2
	float valmax;

	for (i = 0; i < 3; i++) //corr的一列
	{
		float val_max = 0;
		for (j = 0; j < 5696; j++)
		{
			if (fabs(corr[j][i]) > val_max)
			{
				val_max = fabs(corr[j][i]);
				idx[i] = j; //get idx[0], idx[1], idx[2]
			}
		}

		val[i] = val_max; //get val[0], val[1], val[2]
	}

	valmax = 0;
	for (int k = 0; k < 3; k++)
	{
		if (val[k] > valmax)
		{
			valmax = val[k];
			Nid_2 = k;
		}
	}

	std::cout << "val[0] = " << val[0] << endl;
	std::cout << "val[1] = " << val[1] << endl;
	std::cout << "val[2] = " << val[2] << endl;

	std::cout << "idx[0] = " << idx[0] << endl;
	std::cout << "idx[1] = " << idx[1] << endl;
	std::cout << "idx[2] = " << idx[2] << endl;

	if (idx[Nid_2] < 128 * 4)
	{
		PSS_Peak_idx = idx[Nid_2] + Nums_per_subframe * 5;
	}
	else
	{
		PSS_Peak_idx = idx[Nid_2] + 1;
	}

	return Nid_2;
}

//get Nid_1, sf, return Nid_1
int findSSS(float *sss_preprocess_I, float *sss_preprocess_Q, int Nid_2, int pss_peak_idx)
{
	int i, j;
	int Nid_1; //return Nid_1

	int Len_FFT = 128;
	int Len_SSS = 168;
	int sss_peak_idx = pss_peak_idx - Len_FFT * 3 - 10 - 9 * 2 - 1;

	float sum0, sum5;
	float sum0_tmpI, sum0_tmpQ;
	float sum5_tmpI, sum5_tmpQ;
	float sss0_max, sss5_max;
	int sss0_idx, sss5_idx;

	float corr_sss[Len_SSS][2]; //not initialize

	float pss_fI[62];
	float pss_fQ[62];

	float sss_fI[62];
	float sss_fQ[62];

	float pss_tI[Len_FFT]; //not initialize, pss_t[128]
	float pss_tQ[Len_FFT];

	float sss_tI[Len_FFT]; //not initialize
	float sss_tQ[Len_FFT];

	float pss_f_tmpI[Len_FFT]; //not initialize, pss_f_tmp[128]
	float pss_f_tmpQ[Len_FFT];

	float sss_f_tmpI[Len_FFT];
	float sss_f_tmpQ[Len_FFT];

	float ce_I[62];
	float ce_Q[62];

	float tmpI = 0;
	float tmpQ = 0;

	float sss_c_I[62];
	float sss_c_Q[62];

	for (i = 0; i < 128; i++) // get pss_tI[128], pss_tQ[128]
	{
		pss_tI[i] = sss_preprocess_I[pss_peak_idx + i - 1];
		pss_tQ[i] = sss_preprocess_Q[pss_peak_idx + i - 1];
	}

	fft(pss_tI, pss_tQ, Len_FFT, pss_f_tmpI, pss_f_tmpQ);

	//test pss_f_tmp, true
	// for (i = 0; i < 128; i++)
	// {
	// 	cout << pss_f_tmpI[i] << "+" <<pss_f_tmpQ[i] << "j" <<endl;
	// }

	for (i = 0; i < 31; i++)
	{
		pss_fI[i] = pss_f_tmpI[127 - 30 + i]; //pss_fI[0]= pss_f_tmp[127-30]
		pss_fQ[i] = pss_f_tmpQ[127 - 30 + i];
	}

	for (i = 31, j = 1; i < 62; i++, j++)
	{
		pss_fI[i] = pss_f_tmpI[j]; //pss_fI[31] = pss_f_tmpI[1]
		pss_fQ[i] = pss_f_tmpQ[j];
	}

	//test pss_f, true
	// for (i = 0; i < 62; i++)
	// {
	// 	cout << pss_fI[i] << "+" <<pss_fQ[i] << "j" <<endl;
	// }

	genPSS(Nid_2); //get PSSfreqI[62], PSSfreqQ[62]

	//for (i = 0; i < 62; i++)
	//{
	//	cout << PSSfreq_I[i] << "+" <<PSSfreq_Q[i] << "j" <<endl;
	//}

	// //conj(PSSfreq)
	// for (i = 0; i < 62; i++)
	// {
	// 	PSSfreq_Q[i] = 1.0 * PSSfreq_Q[i];
	// }

	for (i = 0; i < 62; i++) //ce_I[62], ce_Q[62]
	{
		//(a + bi)* (c + di) = (a*c - b*d) +(a*d + b*c)i
		ce_I[i] = pss_fI[i] * PSSfreq_I[i] - pss_fQ[i] * PSSfreq_Q[i];
		ce_Q[i] = pss_fI[i] * PSSfreq_Q[i] + pss_fQ[i] * PSSfreq_I[i];
	}

	// test ce, true
	// for (i = 0; i < 62; i++)
	// {
	// 	cout << ce_I[i] << "+" <<ce_Q[i] << "j" <<endl;
	// }

	//sss
	for (i = 0; i < 128; i++) //test sss_I,true
	{
		sss_tI[i] = sss_preprocess_I[sss_peak_idx + i];
		sss_tQ[i] = sss_preprocess_Q[sss_peak_idx + i];
	}

	fft(sss_tI, sss_tQ, Len_FFT, sss_f_tmpI, sss_f_tmpQ);

	for (i = 0; i < 31; i++)  //test sss_f, true
	{
		sss_fI[i] = sss_f_tmpI[127 - 30 + i]; //sss_fI[0] = sss_f_tmpI[127-30]
		sss_fQ[i] = sss_f_tmpQ[127 - 30 + i];
	}

	for (i = 31, j = 1; i < 62; i++, j++)
	{
		sss_fI[i] = sss_f_tmpI[j]; //sss_fI[31] = sss_f_tmpI[1]
		sss_fQ[i] = sss_f_tmpQ[j];
	}

	//conj(ce)
	for (i = 0; i < 62; i++)
	{
		ce_Q[i] = -1.0 * ce_Q[i];
	}

	for (i = 0; i < 62; i++)  // test, true
	{
		//(a + bi)* (c + di) = (a*c - b*d) +( a*d + b*c)i
		sss_c_I[i] = sss_fI[i] * ce_I[i] - sss_fQ[i] * ce_Q[i];
		sss_c_Q[i] = sss_fI[i] * ce_Q[i] + sss_fQ[i] * ce_I[i];
	}
	// test
	// for (i = 0; i < 62; i++)
	// {
	// 	cout << sss_c_I[i] << "+" <<sss_c_Q[i] << "j" <<endl;
	// }

	for (Nid_1 = 0; Nid_1 < 168; Nid_1++)
	{
		genSSS(Nid_1, Nid_2); // SSS_du_0[62], SSS_du_5[62]

		// sss0
		sum0 = 0;
		sum5 = 0;

		sum0_tmpI = 0, sum0_tmpQ = 0;
		sum5_tmpI = 0, sum5_tmpQ = 0;

		for (i = 0; i < 62; i++)
		{
			sum0_tmpI = sum0_tmpI + SSS_du_0[i] * sss_c_I[i];
			sum0_tmpQ = sum0_tmpQ + SSS_du_0[i] * sss_c_Q[i];

			sum5_tmpI = sum5_tmpI + SSS_du_5[i] * sss_c_I[i];
			sum5_tmpQ = sum5_tmpQ + SSS_du_5[i] * sss_c_Q[i];
		}

		sum0 = sqrt(sum0_tmpI * sum0_tmpI + sum0_tmpQ * sum0_tmpQ);
		sum5 = sqrt(sum5_tmpI * sum5_tmpI + sum5_tmpQ * sum5_tmpQ);

		corr_sss[Nid_1][0] = sum0;
		corr_sss[Nid_1][1] = sum5;
	}

	//test corr_sss, true
	// for (int i = 0; i < 168; i++)
	// {
	// 	cout << corr_sss[i][0] << " " << corr_sss[i][1] << endl;
	// }

	//get sss0_max, sss5_max
	sss0_max = 0, sss5_max = 0;
	sss0_idx = 0, sss5_idx = 0;

	for (i = 0; i < 168; i++)
	{
		if (corr_sss[i][0] > sss0_max)
		{
			sss0_max = corr_sss[i][0];
			sss0_idx = i;
		}

		if (corr_sss[i][1] > sss5_max)
		{
			sss5_max = corr_sss[i][1];
			sss5_idx = i;
		}
	}

	if (sss0_max > sss5_max)
	{
		sf = 0;
		Nid_1 = sss0_idx - 1;
	}
	else
	{
		sf = 5;
		Nid_1 = sss5_idx - 1;
	}

	std::cout << "sss0_max = " << sss0_max << " " << "sss5_max = "<< sss5_max <<endl;
	std::cout << "sss0_idx = " << sss0_idx << " " << "sss5_idx = "<< sss5_idx <<endl;

	Nid_1 = Nid_1 + 1;
	return Nid_1;
}

//cell_id = Nid_1*3 + Nid_2
int Set_CellID(int Nid_1, int Nid_2)
{
	return (Nid_1 * 3 + Nid_2);
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

int main()
{
	int count = 0;
	int Nid_2 = 0;
	int Nid_1 = 0;
	int Cell_ID = 0; //Cell_ID = Nid_1 * 3 + Nid_2

	LTESync_init();

	Nid_2 = findPSS(pss_preprocess_I, pss_preprocess_Q); //corr[5696][3]
	Nid_1 = findSSS(sss_preprocess_I, sss_preprocess_Q, Nid_2, PSS_Peak_idx * 2);
	Cell_ID = Set_CellID(Nid_1, Nid_2);

	cout << "Nid_1: " << Nid_1 << endl;
	cout << "Nid_2: " << Nid_2 << endl;
	cout << "Cell_ID: " << Cell_ID << endl;
	cout << "subframe idx: " << sf <<endl;
	cout << "peak_idx: " << PSS_Peak_idx * 2 << endl;

	return 0;
}