#include "PreProcess.h"
#include "Head.h"

struct Signal_short
{
    short I_Data; //2 bytes
    short Q_Data;
};


void PreProcess_init(ifstream &infile)
{
    int CurrentNum = 0;
    GetData(infile); // get I_total[NumSamples_total], Q_total[NumSamples_total]

    while (CurrentNum < NumSamples_total)
    {
        CurrentNum = CurrentNum + NumSamples_subframe;
        Resample(ss_cs_I, ss_cs_Q, ss_cs_I_Dec, ss_cs_Q_Dec);

        conv_same(ss_cs_I_Dec, NumSamples_subframe / 16, filter_ds, filter_ds_length, ss_cs_I_Dec);
        conv_same(ss_cs_Q_Dec, NumSamples_subframe / 16, filter_ds, filter_ds_length, ss_cs_Q_Dec);

        Decimation(ss_cs_I_Dec, NumSamples_subframe / 16, 2, ss_cs_I_process);
        Decimation(ss_cs_Q_Dec, NumSamples_subframe / 16, 2, ss_cs_Q_process);

        ofstream outfile;
        outfile.open("usrp_preprocess.dat", ios::app | ios::binary);

        if (outfile.is_open())
        {
            for (int i = 0; i < CurrentNum / 32; i++)
            {
                outfile.write((char *)(ss_cs_I_process + i), sizeof(float));
                outfile.write((char *)(ss_cs_Q_process + i), sizeof(float));
            }
        }
        outfile.close();
    }
}

void GetData(ifstream &infile) //&infile
{
    Signal_short signal;
    int readBytes = 0;

    if (!infile)
    {
        cout << "Error, can not open the file \n"
             << endl;
        return;
    }

    while (infile.read((char *)&signal, sizeof(signal)) && readBytes < NumSamples_total * 4)
    {
        readBytes += infile.gcount(); //count the read bytes

        if (readBytes % 4 == 0 && readBytes != 0)
        {
            I_total[readBytes / 4 - 1] = signal.I_Data;
            Q_total[readBytes / 4 - 1] = signal.Q_Data;
        }
    }
    //cout << "the data numbers of read is :" << readBytes / 4 << endl;
    
    infile.close();
}

void Resample(float *ss_cs_I, float *ss_cs_Q, float *ss_cs_I_dec, float *ss_cs_Q_dec)
{
    float *ss_cs_I_2Dec = new float[NumSamples_total / 2]; //after the HBF
    float *ss_cs_Q_2Dec = new float[NumSamples_total / 2];

    float *ss_cs_I_4Dec = new float[NumSamples_total / 4]; //after the HBF * 2
    float *ss_cs_Q_4Dec = new float[NumSamples_total / 4];

    float *ss_cs_I_8Dec = new float[NumSamples_total / 8]; //after the HBF * 3
    float *ss_cs_Q_8Dec = new float[NumSamples_total / 8];

    halfband_filter(ss_cs_I, NumSamples_total, hbfiltercoef2, hbfiltercoef2length, ss_cs_I_2Dec);
    halfband_filter(ss_cs_Q, NumSamples_total, hbfiltercoef2, hbfiltercoef2length, ss_cs_Q_2Dec);

    halfband_filter(ss_cs_I_2Dec, NumSamples_total / 2, hbfiltercoef2, hbfiltercoef2length, ss_cs_I_4Dec);
    halfband_filter(ss_cs_Q_2Dec, NumSamples_total / 2, hbfiltercoef2, hbfiltercoef2length, ss_cs_Q_4Dec);

    halfband_filter(ss_cs_I_4Dec, NumSamples_total / 4, hbfiltercoef1, hbfiltercoef1length, ss_cs_I_8Dec);
    halfband_filter(ss_cs_Q_4Dec, NumSamples_total / 4, hbfiltercoef1, hbfiltercoef1length, ss_cs_Q_8Dec);

    halfband_filter(ss_cs_I_8Dec, NumSamples_total / 8, decfiltercoef, decfilterlength, ss_cs_I_Dec);   //after HBF * 4
    halfband_filter(ss_cs_Q_8Dec, NumSamples_total / 8, decfiltercoef, decfilterlength, ss_cs_Q_Dec);

    delete[] ss_cs_I_2Dec;
    delete[] ss_cs_Q_2Dec;

    delete[] ss_cs_I_4Dec;
    delete[] ss_cs_Q_4Dec;

    delete[] ss_cs_I_8Dec;
    delete[] ss_cs_Q_8Dec;
}

void halfband_filter(float *input, int inlength, const float *filtercoef, int filength, float *output)
{
    float *output_temp = new float[inlength];

    conv_same(input, inlength, filter_ds, filter_ds_length, output_temp);
    Decimation(output_temp, inlength, 2, output);

    delete[] output_temp;
}

void Decimation(float *input, int inlength, int Deci, float *output)
{
    int i = 0, k = 0;
    for (i = 0, k = 0; k < inlength; k += Deci, i++)
    {
        output[i] = input[k];
    }
}

void conv(float *x, int x_length, const float *h, int h_length, float *output)
{
    int k = 0, i = 0;
    int Len = x_length + h_length - 1;
    for (k = 0; k < Len; k++)
    {
        for (i = max(0, k + 1 - h_length); i <= min(k, x_length - 1); i++)
        {
            output[k] += x[i] * h[k - i];
        }
    }
}

void conv_same(float *x, int x_length, const float *h, const int h_length, float *output)
{
    int k = 0, i = 0;
    int Len = x_length + h_length - 1;

    float *output_temp = new float[Len]; // output[len]
    return_zero(output_temp, Len);
    // memset(output_temp, 0, Len); 

    for (k = 0; k < Len; k++)
    {
        for (i = max(0, k + 1 - h_length); i <= min(k, x_length - 1); i++)
        {  output_temp[k] += x[i] * h[k - i]; }
    }

    if ((Len - x_length) % 2 == 0)
    {
        for (i = 0, k = (Len - x_length) / 2; i < x_length; i++, k++)
        {
            output[i] = output_temp[k];
        }
    }
    else
    {
        for (i = 0, k = (Len - x_length) / 2 + 1; i < x_length; i++, k++)
        {
            output[i] = output_temp[k];
        }
    }
}

int main(int argc, const char **argv)
{
    ifstream infile("usrp_samples.dat", ios::in | ios::binary);
    PreProcess_init(infile); // usrp_preprocess.dat

    return 0;
}