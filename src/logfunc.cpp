/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "logfunc.h"

// function to generate a vector whose value is equal to its index
// this is useful when plotting a vector because we use the index on the X axis
void genIndexVector(std::vector<float> &x, const int size) {
	x.resize(size, static_cast<float>(0));
	for (int i=0; i<size; i++) {
		x[i] = static_cast<float>(i);
	}
}

// function to be used for logging a float vector in a .dat file (for .gnuplot)
// can be reused for different types of vectors with 32-bit floating point vals
void logVector(const std::string filename, \
	const std::vector<float> &x, \
	const std::vector<float> &y)
{
	// write data in text format to be parsed by gnuplot (change as needed)
	const std::string dat_filename = "../data/" + filename + ".dat";
	std::fstream fd;
	fd.open(dat_filename, std::ios::out);
	fd << "#\tx_axis\ty_axis\n";

	for (unsigned int i = 0; i < x.size(); i++) {
		fd << "\t " << x[i] << "\t";
		// if the number of values on the Y axis is less than on the X tx_axis
		// then we just do not write anything on the Y axis
		if (i < y.size())
			fd << y[i];
		fd << "\n";
	}
	std::cout << "Generated " << dat_filename << " to be used by gnuplot\n";
	fd.close();
}
