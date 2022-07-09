#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

"""
The command-line instructions for recording RF data are only for those
who have the RF dongle (Nooelec NESDR Smart v4 Bundle) available at home.
After you have installed the drivers to work with the RF dongle,
the 8-bit unsigned values for the I/Q pairs can be recorded as follows:

rtl_sdr -f 99.9M -s 2.4M - > iq_samples.raw

The above assumes that we are tuned to the FM station at 99.9 MHz,
we use an RF sample rate of 2.4 Msamples/sec and our file is called
iq_samples.raw (change as you see fit).

For the above use case, the data acquisition runs indefinitely,
hence the recording needs to be stopped by pressing Ctrl+C.
If we wish to stop it after a pre-defined number of samples,
e.g., 12 million I/Q pairs (5 seconds at 2.4 Msamples/sec),
we can use an extra argument:

rtl_sdr -f 99.9M -s 2.4M -n 12000000 - > iq_samples.raw

To check if the raw I/Q data has been recorded properly, place the file
in the "data" sub-folder from your project repo and run this Python file
from the "model" sub-folder. It should produce both the .png image files
(of the PSD estimates) for a couple of blocks, as well as the .wav file.

In the source code below (check lines 90+) you can observe where the
raw_data is read and the normalization of the 8-bit unsigned I/Q samples
to 32-bit float samples (in the range -1 to +1) is done; while the
32-bit floats and the range -1 to +1 are optional choices (used by
many third-party SDR software implementations), it is at the discretion
of each project group to decide how to handle the 8-bit unsigned I/Q samples
in their Python model and C++ implementation.
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import wavfile
from scipy import signal
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD
# for take-home add your functions

# the radio-frequency (RF) sampling rate
# this sampling rate is either configured on RF hardware
# or documented when a raw file with IQ samples is provided
rf_Fs = 2.4e6

# the cutoff frequency to extract the FM channel from raw IQ data
rf_Fc = 100e3

# the number of taps for the low-pass filter to extract the FM channel
# this default value for the width of the impulse response should be changed
# depending on some target objectives, like the width of the transition band
# and/or the minimum expected attenuation from the pass to the stop band
rf_taps = 151

# the decimation rate when reducing the front end sampling rate (i.e., RF)
# to a smaller samping rate at the intermediate frequency (IF) where
# the demodulated data will be split into the mono/stereo/radio data channels
rf_decim = 10

# audio sampling rate (we assume audio will be at 48 KSamples/sec)
audio_Fs = 48e3
# should be the same as rf_Fs / rf_decim / audio_decim

# complete your own settings for the mono channel
# (cutoff freq, audio taps, decimation rate, ...)
audio_Fc = 16e3 #... change as needed (see spec in lab document)
audio_decim = 5 #... change as needed (see spec in lab document)
audio_taps = 151#... change as you see fit

def myLowPassFilter(Fc, Fs, Ntaps):
	h = np.zeros(Ntaps)
	norm_f = Fc/(Fs/2)
	for i in range(Ntaps):
		if (i == (Ntaps -1)/2):
			h[i] = norm_f
		else:
			h[i] = norm_f*np.sin(np.pi*norm_f*(i-(Ntaps-1)/2))/(np.pi*norm_f*(i-(Ntaps-1)/2))
		h[i] = h[i] *(np.sin(i*np.pi/Ntaps))**2
	return h

def convolveFIR(h, input):
	output = np.zeros(len(input) + len(h) -1)
	print(len(input))
	for i in range(len(output)):
		for k in range(len(h)):
			if(i-k>=0 and i-k<len(input)):
				output[i] += h[k]*input[i-k]
	return output


# flag that keeps track if your code is running for
# in-lab (il_vs_th = 0) vs takehome (il_vs_th = 1)
il_vs_th = 1

if __name__ == "__main__":

	# read the raw IQ data from the recorded file
	# IQ data is assumed to be in 8-bits unsigned (and interleaved)
	in_fname = "../data/samples0.raw"
	raw_data = np.fromfile(in_fname, dtype='uint8')
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
	# IQ data is normalized between -1 and +1 in 32-bit float format
	iq_data = (np.float32(raw_data) - 128.0)/128.0
	print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

	# coefficients for the front-end low-pass filter
	rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

	# filter to extract the FM channel (I samples are even, Q samples are odd)
	i_filt = signal.lfilter(rf_coeff, 1.0, iq_data[0::2])
	q_filt = signal.lfilter(rf_coeff, 1.0, iq_data[1::2])

	# downsample the FM channel
	i_ds = i_filt[::rf_decim]
	q_ds = q_filt[::rf_decim]

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# FM demodulator (check the library)
	fm_demod, dummy = fmDemodArctan(i_ds, q_ds)
	# we use a dummy because there is no state for this single-pass model

	# PSD after FM demodulation
	fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0],\
	 		'Demodulated FM (full recording)')

	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle  as for rf_coeff (but different arguments, of course)
		audio_coeff = signal.firwin(audio_taps, audio_Fc/(2.4e5/2), window = 'hann')
		#audio_coeff = np.array([])
	else:
		# to be updated by you for the takehome exercise
		# with your own code for impulse response generation
		audio_coeff = myLowPassFilter(audio_Fc,2.4e5,audio_taps)

	# extract the mono audio data through filtering
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on lfilter
		# same principle as for i_filt or q_filt (but different arguments)
		audio_filt = signal.lfilter(audio_coeff, 1.0, fm_demod)
	else:
		# to be updated by you for the takehome exercise
		# with your own code for single pass convolution
		audio_filt = convolveFIR(audio_coeff, fm_demod)

	# you should uncomment the plots below once you have processed the data

	# PSD after extracting mono audio
	fmPlotPSD(ax1, audio_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Extracted Mono')

	# downsample audio data (see the principle for i_ds or q_ds)
	audio_data = audio_filt[::audio_decim] # to be updated by you during in-lab (same code for takehome)

	# PSD after decimating mono audio
	fmPlotPSD(ax2, audio_data, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')

	# save PSD plots
	fig.savefig("../data/fmMonoBasic.png")
	plt.show()

	# write audio data to file (assumes audio_data samples are -1 to +1)
	out_fname = "../data/fmMonoBasic.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
	# during FM transmission audio samples in the mono channel will contain
	# the sum of the left and right audio channels; hence, we first
	# divide by two the audio sample value and then we rescale to fit
	# in the range offered by 16-bit signed int representation
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")
