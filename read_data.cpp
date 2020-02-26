#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

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

void Get_Short_Data(ifstream &infile, Signal_Short signal);
void Get_Float_Data(ifstream &infile, Signal_Float signal);

int main()
{
    int i;
    Signal_Short short_signal;
    Signal_Float float_signal;

    std::cout << "choose the file to display: " << endl;
    std::cout << "1-------usrp_process_pss.dat" << endl;     //Signal_float
    std::cout << "2-------usrp_process_sss.dat" << endl;     //Signal_float
    std::cout << "3-------usrp_process_PCI_56.dat" << endl;  //Signal_short
    std::cout << "4-------usrp_process_PCI_463.dat" << endl; //Signal_float
    std::cout << "5-------usrp_samples_PCI_147.dat" << endl; //Signal_short

    cin >> i;

    if (i < 1 || i > 5)
    {
        cout << "Error, invalid i:" << i << endl;
        return 0;
    }

    else if (i == 1)
    {
        ifstream infile("/home/leo/vscode/LTECell_Search_V1/usrp_process_pss.dat", ios::in | ios::binary);
        Get_Float_Data(infile, float_signal);
    }

    else if (i == 2)
    {
        ifstream infile("/home/leo/vscode/LTECell_Search_V1/usrp_process_sss.dat", ios::in | ios::binary);
        Get_Float_Data(infile,float_signal);
    }

    else if (i == 3)
    {
        ifstream infile("/home/leo/vscode/LTECell_Search_V1/usrp_samples_PCI_56.dat", ios::in | ios::binary);
        Get_Short_Data(infile, short_signal);
    }

    else if(i == 4)
    {
        ifstream infile("/home/leo/vscode/LTECell_Search_V1/usrp_samples_PCI_463.dat", ios::in | ios::binary);
        Get_Float_Data(infile, float_signal);
    }
    else //if(i == 4)
    {
        ifstream infile("/home/leo/vscode/LTECell_Search_V1/usrp_samples_PCI_147.dat", ios::in | ios::binary);
        Get_Short_Data(infile, short_signal);
    }

    //cout << "the numbers of total data is: " << i << endl;
    return 0;
}

void Get_Short_Data(ifstream &infile, Signal_Short signal) //&infile
{
    int readBytes = 0;
    int count = 0;

    if (!infile)
    {
        cout << "Error, can not open the file \n"
             << endl;
        return;
    }

    while (infile.read((char *)&signal, sizeof(signal)))
    {
        readBytes += infile.gcount(); //count the read bytes
        cout << signal.I_Data << " " << signal.Q_Data << endl;
        count++;
    }

    cout << "the bytes numbers of read is : " << readBytes << endl;
    cout << "the  data numbers of read is : " << count << endl;

    infile.close();
}


void Get_Float_Data(ifstream &infile, Signal_Float signal) //&infile
{
    int readBytes = 0;
    int count = 0;

    if (!infile)
    {
        cout << "Error, can not open the file \n"
             << endl;
        return;
    }

    while (infile.read((char *)&signal, sizeof(signal)))
    {
        readBytes += infile.gcount(); //count the read bytes
        cout << signal.I_Data << " " << signal.Q_Data << endl;
        count++;
    }

    cout << "the bytes numbers of read is : " << readBytes << endl;
    cout << "the  data numbers of read is : " << count << endl;

    infile.close();
}