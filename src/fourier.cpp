/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

// source code for Fourier-family of functions
#include "dy4.h"
#include "fourier.h"

// just DFT function (no FFT yet)
void DFT(const std::vector<float> &x, std::vector<std::complex<float>> &Xf) {
	Xf.resize(x.size(), static_cast<std::complex<float>>(0));
	for (unsigned int m = 0; m < Xf.size(); m++) {
		for (unsigned int k = 0; k < x.size(); k++) {
				std::complex<float> expval(0, -2*PI*(k*m) / x.size());
				Xf[m] += x[k] * std::exp(expval);
		}
	}
}

// function to compute the magnitude values in a complex vector
void computeVectorMagnitude(const std::vector<std::complex<float>> &Xf, std::vector<float> &Xmag)
{
	// only the positive frequencies
	Xmag.resize(Xf.size(), static_cast<float>(0));
 	for (unsigned int i = 0; i < Xf.size(); i++) {
		Xmag[i] = std::abs(Xf[i])/Xf.size();
	}
}

// add your own code to estimate the PSD
// NFFT is already definded in dy4.h
void estimatePSD(const std::vector<float> &samples, const float Fs, std::vector<float> &freq, std::vector<float> &psd_est) {

	float df = Fs/NFFT;

	//initialize frequency
	freq.clear(); freq.resize(Fs/2/df,  static_cast<float>(0));
	for(unsigned int i = 0; i<freq.size();i++)
		freq[i] = i*df;

	// design the Hann window used to smoothen the discrete data in order
	// to reduce the spectral leakage after the Fourier transform
	std::vector<float> hann;
	hann.resize(NFFT, static_cast<float>(0));
	for (unsigned int i = 0; i<NFFT; i++)
		hann[i] = std::pow(std::sin(i*PI/NFFT),2);
	// create an empty list where the PSD for each segment is computed
	std::vector<float> psd_list;

	// samples should be a multiple of frequency bins, so
	// the number of segments used for estimation is an integer
	// note: for this to work you must provide an argument for the
	// number of frequency bins not greater than the number of samples!
	unsigned int no_segments = (int)std::floor(samples.size()/(float)NFFT);
	std::vector<std::complex<float>>  Xf;
	std::vector<float> psd_seg, windowed_samples;

	// iterate through all the segments
	for (unsigned int k = 0; k<no_segments ; k++){
		windowed_samples.clear();windowed_samples.resize((int)NFFT, static_cast<float>(0));

		// apply the hann window (using pointwise multiplication)
		// before computing the Fourier transform on a segment
		for (unsigned int i = 0; i<NFFT;i++){
			windowed_samples[i] = samples[(int)k*NFFT+i] * hann[i];
		}

		// compute the Fourier transform using the built-in FFT from numpy
		Xf.clear();
		DFT(windowed_samples, Xf); //check if needed FFT_improved

		// since input is real, we keep only the positive half of the spectrum
		// however, we will also add the signal energy of negative frequencies
		// to have a better a more accurate PSD estimate when plotting
		std::vector<std::complex<float>>::iterator end = Xf.begin()+(int)NFFT/2;
		Xf.assign(Xf.begin(), end); // keep only positive freq bins

		//psd_seg
		psd_seg.clear(); psd_seg.resize(Xf.size(),  static_cast<float>(0));
		for (unsigned int i = 0; i<Xf.size();i++){
			psd_seg[i] = 2*(1/(Fs*NFFT/2)) * std::pow(std::abs(Xf[i]),2) ;// compute signal power
			 // add the energy from the negative freq bins
			 // translate to the decibel (dB) scale
			 psd_seg[i] = 10*std::log10(psd_seg[i]);
		}

		// append to the list where PSD for each segment is stored
		// in sequential order (first segment, followed by the second one, ...)
		psd_list.insert(psd_list.end(), psd_seg.begin(), psd_seg.end());

	}
	//printRealVector1(psd_list);

	// compute the estimate to be returned by the function through averaging
	psd_est.clear();
	psd_est.resize((int)NFFT/2,  static_cast<float>(0));

	// iterate through all the frequency bins (positive freq only)
	// from all segments and average them (one bin at a time ...)
	for (unsigned int k = 0; k< (int)(NFFT/2); k++){
		// iterate through all the segments
		for( unsigned int l = 0; l< no_segments; l++){
			psd_est[k] += psd_list[k + l*(int)(NFFT/2)];
		}
		// compute the estimate for each bin
		psd_est[k] = psd_est[k] / (float)no_segments;
	}
	//printRealVector1(psd_est);

}

// added IDFT

void IDFT(const std::vector<std::complex<float>> &Xf, std::vector<std::complex<float>> &x) {
	x.resize(Xf.size(), static_cast<std::complex<float>>(0));
	for (unsigned int k = 0; k < x.size(); k++) {
		for (unsigned int m = 0; m < x.size(); m++) {
			std::complex<float> expval(0, 2*PI*(k*m) / Xf.size());
			x[k] += Xf[m] * std::exp(expval);
		}
		x[k] /= Xf.size();
	}
}

// added FFT

unsigned int swap_bits(unsigned int x, unsigned char i, unsigned char j) {

  unsigned char bit_i = (x >> i) & 0x1L;
  unsigned char bit_j = (x >> j) & 0x1L;

  unsigned int val = x;
  val = bit_i ? (val | (0x1L << j)) : (val & ~(0x1L << j));
  val = bit_j ? (val | (0x1L << i)) : (val & ~(0x1L << i));

  return val;
}

unsigned int bit_reversal(unsigned int x, unsigned char bit_size) {

  unsigned int val = x;

  for (int i=0; i < int(bit_size/2); i++)
    val = swap_bits(val, i, bit_size-1-i);

  return val;
}

void compute_twiddles(std::vector<std::complex<float>> &twiddles) {
  for (int k=0; k<(int)twiddles.size(); k++) {
      std::complex<float> expval(0.0, -2*PI*float(k)/ NFFT);
      twiddles[k] = std::exp(expval);
  }
}

void FFT_recursive(const std::vector<std::complex<float>> &x, \
  std::vector<std::complex<float>> &Xf) {

  if (x.size() > 1) {
    // declare vectors and allocate space for the even and odd halves
    std::vector<std::complex<float>> xe(int(x.size()/2)), xo(int(x.size()/2));
    std::vector<std::complex<float>> Xfe(int(x.size()/2)), Xfo(int(x.size()/2));

    // split into even and odd halves
    for (int k=0; k<(int)x.size(); k++)
      if ((k%2) == 0) xe[k/2] = x[k];
      else xo[k/2] = x[k];

    // call recursively FFT of half size for even and odd halves respectively
    FFT_recursive(xe, Xfe);
    FFT_recursive(xo, Xfo);

    // merge the results from the odd/even FFTs (each of half the size)
    for (int k=0; k<(int)xe.size(); k++) {
        std::complex<float> expval(0.0, -2*PI*float(k)/ x.size());
        std::complex<float> twiddle = std::exp(expval);
        Xf[k]           = Xfe[k] + twiddle * Xfo[k];
        Xf[k+xe.size()] = Xfe[k] - twiddle * Xfo[k];
    }
  } else {
    // end of recursion - copy time domain samples to frequency bins (default values)
    Xf[0] = x[0];
  }
}

void FFT_improved(const std::vector<std::complex<float>> &x, \
  std::vector<std::complex<float>> &Xf, \
  const std::vector<std::complex<float>> &twiddles, \
  const unsigned char recursion_level) {

  if (x.size() > 1) {
    int half_size = int(x.size()/2);
    std::vector<std::complex<float>> xe(half_size), xo(half_size);
    std::vector<std::complex<float>> Xfe(half_size), Xfo(half_size);

    for (int k=0; k<half_size; k++) {
      xe[k] = x[k*2];
      xo[k] = x[k*2+1];
    }

    FFT_improved(xe, Xfe, twiddles, recursion_level+1);
    FFT_improved(xo, Xfo, twiddles, recursion_level+1);

    for (int k=0; k<half_size; k++) {
        Xf[k]           = Xfe[k] + twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
        Xf[k+half_size] = Xfe[k] - twiddles[k*(1<<(recursion_level-1))] * Xfo[k];
    }
  } else {
    Xf[0] = x[0];
  }
}

void FFT_optimized(const std::vector<std::complex<float>> &x, \
  std::vector<std::complex<float>> &Xf, \
  const std::vector<std::complex<float>> &twiddles) {

  unsigned char no_levels = (unsigned char)std::log2((float)x.size());
  for (unsigned int i=0; i<x.size(); i++) {
    Xf[i] = x[bit_reversal(i, no_levels)];
  }

  unsigned int step_size = 1;

  std::complex<float> tmp;
  for (unsigned char l=0; l<no_levels; l++) {
    for (unsigned int p=0; p<x.size(); p+=2*step_size) {
      for (unsigned int k=p; k<p+step_size; k++) {
        tmp             = Xf[k] + twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
        Xf[k+step_size] = Xf[k] - twiddles[(k-p)*(1<<(no_levels-1-l))] * Xf[k+step_size];
        Xf[k]           = tmp;
      }
    }
    step_size *= 2;
  }
}
