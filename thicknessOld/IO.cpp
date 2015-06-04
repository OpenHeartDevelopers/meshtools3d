#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <set>
#include <list>

using namespace std;

// Reads-in a temporal file of integers with suffix .dat with specified length
void readTemporalIntegerDataFile(string file, unsigned int length, int* &data, int t)
{
	string dat = ".dat";
	stringstream ss;
	string str_t;
	ss << t;
	ss >> str_t;
	string total_datafilename = file + str_t + dat;

	ifstream datafile(total_datafilename.c_str());

	if (datafile.is_open())
	{
		data = new int [length];

		for(unsigned int i=0;i<length;i++)
		{
			datafile >> data[i];
			//if(data[i] != 0)
				//cout << data[i] << "\n";
		}
	}
	else
	{
		cerr << "Can't Open " << datafile << "!" << endl << flush;
	}
}

// Reads-in temporal file of integers with specified ending and length given by first line of file
void readTemporalIntegerListFile(string file, string end, unsigned int &length, int* &data, int t)
{
	stringstream ss;
	string str_t;
	ss << t;
	ss >> str_t;
	string total_datafilename = file + str_t + end;

	ifstream datafile(total_datafilename.c_str());

	if (datafile.is_open())
	{
		datafile >> length;
		
		data = new int [length];

		for(unsigned int i=0;i<length;i++)
		{
			datafile >> data[i];
		} 
	}
	else
	{
		cerr << "Can't Open " << datafile << "!" << endl << flush;
	}
}

// Reads-in a temporal file of doubles with the suffix .dat and specified length
void readTemporalDataFile(const char* file, unsigned int length, double* &data, int t)
{
	string dat = ".dat";
	stringstream ss;
	string str_t;
	ss << t;
	ss >> str_t;
	string total_datafilename = file + str_t + dat;
	ifstream datafile(total_datafilename.c_str());

	if (datafile.is_open())
	{
		data = new double [length];

		for(unsigned int i=0;i<length;i++)
		{
			datafile >> data[i];
		}
	}
	else
	{
		cerr << "Can't Open " << datafile << "!" << endl << flush;
	}
	datafile.close();
}

// Reads-in a data file of doubles with the suffix .dat and specified length
void readDataFile(const char* file, unsigned int length, double* &data)
{
	ifstream datafile(file);
	
	if (datafile.is_open())
	{
		data = new double [length];

		for(unsigned int i=0;i<length;i++)
		{
			datafile >> data[i];
		}
	}
	else
	{
		cerr << "Can't Open " << datafile << "!" << endl << flush;
	}
	datafile.close();
}

// Reads-in an integer data file of a specified length
void readIntegerDataFile(const char* file, unsigned int length, int* &data)
{
        ifstream datafile(file);

        if (datafile.is_open())
        {
                data = new int [length];

                for(unsigned int i=0;i<length;i++)
                {
                        datafile >> data[i];
                }
        }
        else
        {
                cerr << "Can't Open " << datafile << "!" << endl << flush;
        }
        datafile.close();
}


// Writes-out a temporal file of integers with a specified ending with the first line giving the file length
void writeTemporalIntegerListFile(string file, string end, unsigned int length, int* data, int t)
{
	stringstream ss;
	string str_t;
	ss << t;
	ss >> str_t;
	string total_datafilename = file + str_t + end;

	ofstream datafile(total_datafilename.c_str());

	if (datafile.is_open())
	{
		datafile << length << "\n";
		
		for(unsigned int i=0;i<length;i++)
		{
			datafile << data[i] << "\n";
		}
	}
	else
	{
		cerr << "Can't Open " << datafile << "!" << endl << flush;
	}
}

// Writes-out a temporal file of integers with suffix .dat and specified length
void writeTemporalIntegerDataFile(string file, unsigned int length, int* data, int t)
{
	string dat = ".dat";
	stringstream ss;
	string str_t;
	ss << t;
	ss >> str_t;
	string total_datafilename = file + str_t + dat;
	ofstream datafile(total_datafilename.c_str());

	for(unsigned int i=0;i<length;i++)
	{
		datafile << data[i] << "\n";
	}

	datafile.close();
}

// Writes-out temporal file of doubles with suffix .dat and of specified length
void writeTemporalDataFile(string file, unsigned int length, double* data, int t)
{
	string dat = ".dat";
	stringstream ss;
	string str_t;
	ss << t;
	ss >> str_t;
	string total_datafilename = file + str_t + dat;
	ofstream datafile(total_datafilename.c_str());
	
		for(unsigned int i=0;i<length;i++)
		{
			datafile << data[i] << "\n";
		}
		
		datafile.close();
}

// Writes-out an integer file of specified length
void writeIntegerDataFile(string file, unsigned int length, int* data)
{
	ofstream datafile(file.c_str());

	for(unsigned int i=0;i<length;i++)
	{		
		datafile << data[i] << "\n";
	}
		
		datafile.close();
}

// Writes-out an integer data file of specified length
 void writeIntegerListFile(string file, unsigned int length, int* data)
{
	ofstream datafile(file.c_str());
	
	datafile << length << "\n";

        for(unsigned int i=0;i<length;i++)
        {
	         datafile << data[i] << "\n";
        }
        datafile.close();
}

// Writes-out a data file of specified length
void writeDataFile(string file, unsigned int length, double* data)
{
	ofstream datafile(file.c_str());
	
	for(unsigned int i=0;i<length;i++)
	{		
		datafile << data[i] << "\n";
	}
		
		datafile.close();
}

