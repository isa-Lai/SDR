/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_FILTER_H
#define DY4_FILTER_H

// add headers as needed
#include <iostream>
#include <vector>
#include <cmath>

void fmRRC(float, unsigned short int, std::vector<float> &);
//pllIn
void fmPll(float & , float &, float & , float &, float &,float &, std::vector<float> &, std::vector<float> &,  float , float , float ncoScale = 1.0, float phaseAdjust = 0.0, float normBandwidth = 0.01);
void RDS_fmPll(float &, float &, float &,float & , float &, float & , float &,  std::vector<float> &,std::vector<float> &, std::vector<float> &,  float , float , float ncoScale = 1.0, float phaseAdjust = 0.0, float normBandwidth = 0.001);
// declaration of a function prototypes
void allpass(std::vector<float> &, std::vector<float> &, std::vector<float> &);
void impulseResponseLPF(float, float, unsigned short int, std::vector<float> &);
void impulseResponseBPF(float, float, float, unsigned short int, std::vector<float> &);
void impulseResponseLPF_amp(float, float, unsigned short int, int, std::vector<float> &);
void convolveFIR(std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void convolution_with_state(std::vector<float> &, std::vector<float> &, const std::vector<float> &, const std::vector<float> &);
void conv_ds(std::vector<float> &, const std::vector<float> &, std::vector<float> &, const std::vector<float> &,const std::vector<float> &, const int, std::vector<float> &, std::vector<float> &);
//void conv_ds(std::vector<float> &, const std::vector<float> &,const std::vector<float> &, const int, std::vector<float> &);
void resampling(std::vector<float> &, const std::vector<float> &, const std::vector<float> &, const int, const int, std::vector<float> &);



#endif // DY4_FILTER_H
