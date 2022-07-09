/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#ifndef DY4_LOGFUNC_H
#define DY4_LOGFUNC_H

// add headers as needed
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

void genIndexVector(std::vector<float> &, \
  const int);

void logVector(const std::string, \
  const std::vector<float> &, \
  const std::vector<float> &);

#endif // DY4_LOGFUNC_H
