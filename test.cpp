#include <iostream>
#include <math.h>
#include <string.h>
using namespace std;

#define pi 3.14159265
float PSSfreq_I[62];
float PSSfreq_Q[62];

float SSS_du_0[62];
float SSS_du_5[62];

void genPSS(int);
void genSSS(int, int);
int mod(float a, int b);
int mod(int a, int b);

int main()
{
	genPSS(1);
	genSSS(150,1);

	for (int k = 0; k < 62; k++)
	{
		cout << SSS_du_0[k] << " " << SSS_du_5[k] <<endl;
	}

	return 0;
}
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
	m1 = mod(float(m0 + floor(m_prime / 31) + 1), 31);

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