#include "PreProcess.h"
#include <cstdio> //C++

void PreProcess_init()
{
    Signal_Short short_signal;  //read short data
    Signal_Float float_signal;  //read float data

    const char *filePath1 = "/home/leo/vscode/LTECell_Search/usrp_process_pss.dat";
    const char *filePath2 = "/home/leo/vscode/LTECell_Search/usrp_process_sss.dat";

    //ifstream infile("/home/leo/vscode/LTECell_Search/usrp_samples_PCI_463.dat", ios::in | ios::binary);
    //Get_Float_Data(infile); // get I_total, Q_total

    ifstream infile("usrp_samples_PCI_56.dat", ios::in | ios::binary);
    Get_Short_Data(infile); // get I_total, Q_total

    fstream infile1, infile2;
    infile1.open("/home/leo/vscode/LTECell_Search/usrp_process_pss.dat", ios::in);
    infile2.open("/home/leo/vscode/LTECell_Search/usrp_process_sss.dat", ios::in);

    if (!infile1 && !infile2) //if infile1 and infile2 don't exist
    {
        cout << "creat the file" << endl;
        WriteData();
    }

    else // if infile1 exists or infile2 exists
    {
        if (infile1.is_open())
        {
            infile1.close();

            if (remove(filePath1) == 0)
                cout << "delete the original file: usrp_process_pss.dat" << endl;
        }

        if (infile2.is_open())
        {
            infile2.close();

            if (remove(filePath2) == 0)
                cout << "delete the original file: usrp_process_sss.dat" << endl;
        }

        cout << "creat the file" << endl;
        WriteData();
    }
}

void WriteData()
{
    int i = 0;
    int j = 0;
    int CurrentNum = 0; // the numbers data of read now

    while (CurrentNum < NumSamples_total)
    {
        for (i = 0, j = 0; i < NumSamples_subframe; i++, j++)
        {
            ss_cs_I[i] = I_total[j + CurrentNum];
            ss_cs_Q[i] = Q_total[j + CurrentNum];
        }
        CurrentNum = CurrentNum + NumSamples_subframe;
        //Loop 1, CurrentNum:   0   --> 30720
        //Loop 2, CurrentNum: 30720 --> 64140
        //Loop 3, CurrentNum: 64140 --> 92160
        //Loop 4, CurrentNum: 92160 --> 122880
        //Loop 5, CurrentNum:122880 --> 153600
        //Loop 6, CurrentNum:153600 --> 184320

        Resample(ss_cs_I, ss_cs_Q, ss_cs_I_Dec, ss_cs_Q_Dec);

        ofstream outfile1;
        outfile1.open("usrp_process_sss.dat", ios::app | ios::binary);

        if (outfile1.is_open())
        {
            for (int i = 0; i < NumSamples_subframe / 16; i++)
            {
                outfile1.write((char *)(ss_cs_I_Dec + i), sizeof(float));
                outfile1.write((char *)(ss_cs_Q_Dec + i), sizeof(float));
            }
        }

        conv_same(ss_cs_I_Dec, NumSamples_subframe / 16, filter_ds, filter_ds_length, ss_cs_I_Dec);
        conv_same(ss_cs_Q_Dec, NumSamples_subframe / 16, filter_ds, filter_ds_length, ss_cs_Q_Dec);

        Decimation(ss_cs_I_Dec, NumSamples_subframe / 16, 2, ss_cs_I_process);
        Decimation(ss_cs_Q_Dec, NumSamples_subframe / 16, 2, ss_cs_Q_process);

        ofstream outfile2;
        outfile2.open("usrp_process_pss.dat", ios::app | ios::binary);

        if (outfile2.is_open())
        {
            for (int i = 0; i < NumSamples_subframe / 32; i++)
            {
                outfile2.write((char *)(ss_cs_I_process + i), sizeof(float));
                outfile2.write((char *)(ss_cs_Q_process + i), sizeof(float));
            }
        }

        outfile1.close();
        outfile2.close();
    }
}

void Get_Short_Data(ifstream &infile) //&infile
{
    Signal_Short signal;
    int readBytes = 0;

    if (!infile)
    {
        cout << "Error, can not open the file \n"
             << endl;
        return;
    }

    while (infile.read((char *)&signal, sizeof(signal)) && readBytes < NumSamples_total * 4)
    {
        if (readBytes % 4 == 0)
        {
            I_total[readBytes / 4] = float(signal.I_Data) / 1024;
            Q_total[readBytes / 4] = float(signal.Q_Data) / 1024;
        }

        readBytes += infile.gcount(); //count the read bytes
    }
    //cout << "the data numbers of read is :" << readBytes / 4 << endl;

    infile.close();
}

void Get_Float_Data(ifstream &infile) //&infile
{
    Signal_Float signal;
    int readBytes = 0;

    if (!infile)
    {
        cout << "Error, can not open the file \n"
             << endl;
        return;
    }

    while (infile.read((char *)&signal, sizeof(signal)) && readBytes < NumSamples_total * 8)
    {
        if (readBytes % 8 == 0)
        {
            I_total[readBytes / 8] = signal.I_Data / 1024;
            Q_total[readBytes / 8] = signal.Q_Data / 1024;
            //cout << signal.I_Data << " " << signal.Q_Data << endl;
        }

        readBytes += infile.gcount(); //count the read bytes
    }
    //cout << "the data numbers of read is :" << readBytes / 4 << endl;

    infile.close();
}

void Resample(float *ss_cs_I, float *ss_cs_Q, float *ss_cs_I_dec, float *ss_cs_Q_dec)
{
    float *ss_cs_I_2Dec = new float[NumSamples_subframe / 2]; //after the HBF
    float *ss_cs_Q_2Dec = new float[NumSamples_subframe / 2];

    float *ss_cs_I_4Dec = new float[NumSamples_subframe / 4]; //after the HBF * 2
    float *ss_cs_Q_4Dec = new float[NumSamples_subframe / 4];

    float *ss_cs_I_8Dec = new float[NumSamples_subframe / 8]; //after the HBF * 3
    float *ss_cs_Q_8Dec = new float[NumSamples_subframe / 8];

    halfband_filter(ss_cs_I, NumSamples_subframe, hbfiltercoef2, hbfiltercoef2length, ss_cs_I_2Dec);
    halfband_filter(ss_cs_Q, NumSamples_subframe, hbfiltercoef2, hbfiltercoef2length, ss_cs_Q_2Dec);

    halfband_filter(ss_cs_I_2Dec, NumSamples_subframe / 2, hbfiltercoef2, hbfiltercoef2length, ss_cs_I_4Dec);
    halfband_filter(ss_cs_Q_2Dec, NumSamples_subframe / 2, hbfiltercoef2, hbfiltercoef2length, ss_cs_Q_4Dec);

    halfband_filter(ss_cs_I_4Dec, NumSamples_subframe / 4, hbfiltercoef1, hbfiltercoef1length, ss_cs_I_8Dec);
    halfband_filter(ss_cs_Q_4Dec, NumSamples_subframe / 4, hbfiltercoef1, hbfiltercoef1length, ss_cs_Q_8Dec);

    halfband_filter(ss_cs_I_8Dec, NumSamples_subframe / 8, decfiltercoef, decfilterlength, ss_cs_I_Dec); //after HBF * 4
    halfband_filter(ss_cs_Q_8Dec, NumSamples_subframe / 8, decfiltercoef, decfilterlength, ss_cs_Q_Dec);

    delete[] ss_cs_I_2Dec;
    delete[] ss_cs_Q_2Dec;

    delete[] ss_cs_I_4Dec;
    delete[] ss_cs_Q_4Dec;

    delete[] ss_cs_I_8Dec;
    delete[] ss_cs_Q_8Dec;
}

void halfband_filter(float *input, int inlength, const float *filtercoef, int filength, float *output)
{
    float *output_hbf = new float[inlength];

    conv_same(input, inlength, filter_ds, filter_ds_length, output_hbf);
    Decimation(output_hbf, inlength, 2, output);

    delete[] output_hbf;
}

void Decimation(float *input, int inlength, int Deci, float *output)
{
    int i = 0, k = 0;
    for (i = 0, k = 0; k < inlength; k += Deci, i++)
    {
        output[i] = input[k];
    }
}

void conv_same(float *x, int x_length, const float *h, const int h_length, float *output)
{
    int k = 0, i = 0;
    int Len = x_length + h_length - 1;

    float *output_temp = new float[Len]; // output_temp[len]
    //return_zero(output_temp, Len);
    memset(output_temp, 0, Len);

    for (k = 0; k < Len; k++)
    {
        for (i = max(0, k + 1 - h_length); i <= min(k, x_length - 1); i++)
        {
            output_temp[k] += x[i] * h[k - i];
            //output_temp[k] = 2;
        }
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
    delete[] output_temp;
}

int main()
{
    int count = 0; 
    PreProcess_init(); // usrp_preprocess.dat

    //cout << "count = " << count << endl;
    return 0;
}