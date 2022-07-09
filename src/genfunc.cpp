/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "genfunc.h"

// some basic support functions for signal generation from the previous lab
void generateSin(std::vector<float> &t, std::vector<float> &x, float Fs, float interval, float frequency = 7.0, float amplitude = 5.0, float phase = 0.0)
{
	t.resize(0); x.resize(0);
	float dt = 1/Fs;
	for (float i = 0.0; i < interval; i += dt) {
		t.push_back(i);
		x.push_back(amplitude*std::sin(2*PI*frequency*i+phase));
	}
}

void addSin(const std::vector<std::vector<float>> &sv, std::vector<float> &added)
{
	for (unsigned int i = 0.0; i < sv[0].size(); i ++) {
		float addval = 0.0;
		for (unsigned int k = 0; k < sv.size(); k++)
			addval += sv[k][i];
		added.push_back(addval);
	}
}

void generateRandomSamples(std::vector<float> &x, unsigned int N, unsigned short int max, unsigned char precision)
{
	x.resize(N);
	int int_radom_max = 2*(max * static_cast<int>(pow(10,precision)));
	for (unsigned int i = 0; i < x.size(); i++) {
		x[i] = (static_cast<float>(std::rand() % int_radom_max));
		x[i] = (x[i] - (int_radom_max/2))/pow(10,precision);
	}
}
