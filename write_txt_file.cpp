#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;

void GetData(ifstream &infile);

int main()
{
    ifstream infile("/home/leo/vscode/LTECell_Search_V1/usrp_process_sss.dat", ios::in | ios::binary);
    GetData(infile);
    //cout << "the numbers of total data is: " << i << endl;
    return 0;
}

void GetData(ifstream &infile) //&infile
{
    struct Signal_Float
    {
        float I_Data; //2 bytes
        float Q_Data;
    };

    Signal_Float signal;
    int readBytes = 0;
    int count = 0;

    ofstream outfile;
    outfile.open("PCI_404_SSS.txt", ios::out);

    if (!infile)
    {
        cout << "Error, can not open the file \n"
             << endl;
        return;
    }

    // read the signal, then write to the file
    while (infile.read((char *)&signal, sizeof(signal)))
    {
        readBytes += infile.gcount(); //count the read bytes
        outfile << signal.I_Data << " " << signal.Q_Data << endl;
    }

    cout << "the data numbers of read is : " << readBytes / 8 << endl;

    infile.close();
    outfile.close();
}