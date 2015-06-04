#ifndef _IO_H_
#define _IO_H_

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

// Function to load-in temporal data files (e.g. vm.0.dat)
void readTemporalDataFile(const char* file, unsigned int length, double* &data, int t);

// Function to load-in temporal list file (first line gives length)
void readTemporalIntegerListFile(string file, string end, unsigned int &length, int* &data, int t);

// Function to load-in normal single data file
void readDataFile(const char* file, unsigned int length, double* &data);

// Function to load-in normal single data file
void readIntegerDataFile(const char* file, unsigned int length, int* &data);

// Function to write temporal list file (first line gives length)
void writeTemporalIntegerListFile(string file, string end, unsigned int length, int* data, int t);

// Function to load-in temporal integer data files (e.g. vm.0.dat)
void readTemporalIntegerDataFile(string file, unsigned int length, int* &data, int t);

// Function to write-out temporal integer data files (e.g. vm.0.dat)
void writeTemporalIntegerDataFile(string file, unsigned int length, int* data, int t);

// Function to write-out temporal data files (e.g. vm.0.dat)
void writeTemporalDataFile(string file, unsigned int length, double* data, int t);

// Function to write-out data files
void writeIntegerDataFile(string file, unsigned int length, int* data);

// Function to write-out data files
void writeIntegerListFile(string file, unsigned int length, int* data);

// Function to write-out data files
void writeDataFile(string file, unsigned int length, double* data);


#endif
