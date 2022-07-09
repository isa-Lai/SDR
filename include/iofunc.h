/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_IOFUNC_H
#define DY4_IOFUNC_H

// add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <complex>

// declaration of a function prototypes
void printRealVector(const std::vector<float> &);
void printRealVectorInt(const std::vector<int> &);

void printComplexVector(const std::vector<std::complex<float>> &);

void readStdinBlockData(unsigned int, unsigned int, std::vector<float> &);

void readBinData(const std::string, std::vector<float> &);

void readRawFile(const std::string, std::vector<uint8_t> &);

void writeBinData(const std::string, const std::vector<float> &);

#endif // DY4_IOFUNC_H
