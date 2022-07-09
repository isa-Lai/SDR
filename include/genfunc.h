/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_GENFUNC_H
#define DY4_GENFUNC_H

// add headers as needed
#include <iostream>
#include <vector>
#include <complex>

// declaration of a function prototypes
void generateSin(std::vector<float> &, std::vector<float> &, float, \
	float, float, float, float);

void addSin(const std::vector<std::vector<float>> &, std::vector<float> &);

void generateRandomSamples(std::vector<float> &, unsigned int, \
	unsigned short int, unsigned char);

#endif // DY4_GENFUNC_H
