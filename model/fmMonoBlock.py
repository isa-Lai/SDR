#
# Comp Eng 3DY4 (Computer Systems Integration Project)
#
# Copyright by Nicola Nicolici
# Department of Electrical and Computer Engineering
# McMaster University
# Ontario, Canada
#

import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math

# use fmDemodArctan and fmPlotPSD
from fmSupportLib import fmDemodArctan, fmPlotPSD, myDemod
# for take-home add your functions

rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 151
rf_decim = 10

audio_Fs = 48e3
audio_decim = 800
up = 147
# add other settings for audio, like filter taps, ...
audio_taps = 151
audio_Fc = 16e3

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

def block_processing_with_state(h,xb,state):
	yb = np.zeros(len(xb))
	for n in range(len(yb)):
		for k in range(len(h)):
			if n-k>=0:
				yb[n] += h[k] *xb[n-k]
			#to compute with the initial value
			else:
				yb[n] += h[k] * state[n-k]

	if(len(xb)<len(state)):
		# original values are push to the front
		# values in xb is inserted at the back
		state[0:len(state)-len(xb)] = state[len(xb):]
		state[len(state)-len(xb):] = xb
	else:
		# since xb have enough value stored in state
		# direcly assign
		state = xb[len(xb)-len(state):]
	return yb,state

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
	rf_coeff = myLowPassFilter(rf_Fc, rf_Fs, rf_taps)
	

	# coefficients for the filter to extract mono audio
	if il_vs_th == 0:
		# to be updated by you during the in-lab session based on firwin
		# same principle  as for rf_coeff (but different arguments, of course)
		audio_coeff = signal.firwin(audio_taps, audio_Fc/(2.4e5/2), window=('hann'))
	else:
		# to be updated by you for the takehome exercise
		# with your own code for impulse response generation
		audio_coeff = up*myLowPassFilter(audio_Fc,2.4e5*up,audio_taps*up)

	# set up the subfigures for plotting
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)

	# select a block_size that is a multiple of KB
	# and a multiple of decimation factors
	block_size = 1024 * rf_decim * audio_decim * 2
	block_count = 0

	# states needed for continuity in block processing
	state_i_lpf_100k = np.zeros(rf_taps-1)
	state_q_lpf_100k = np.zeros(rf_taps-1)
	state_audio_lpf_100k = np.zeros(audio_taps*up-1)
	state_phase = 0
	state_Q  = 0
	# add state as needed for the mono channel filter

	# audio buffer that stores all the audio blocks
	audio_data = np.array([]) # used to concatenate filtered blocks (audio data)

	# if the number of samples in the last block is less than the block size
	# it is fine to ignore the last few samples from the raw IQ file
	while (block_count+1)*block_size < len(iq_data):

		# if you wish to have shorter runtimes while troubleshooting
		# you can control the above loop exit condition as you see fit
		print('Processing block ' + str(block_count))

		# filter to extract the FM channel (I samples are even, Q samples are odd)
		i_filt, state_i_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size:(block_count+1)*block_size:2],
				zi=state_i_lpf_100k)
		q_filt, state_q_lpf_100k = signal.lfilter(rf_coeff, 1.0, \
				iq_data[(block_count)*block_size+1:(block_count+1)*block_size:2],
				zi=state_q_lpf_100k)

		# downsample the I/Q data from the FM channel
		i_ds = i_filt[::rf_decim]
		q_ds = q_filt[::rf_decim]

		# FM demodulator
		if il_vs_th == 0:
			# already given to you for the in-lab
			# take particular notice of the "special" state-saving
			fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)
		else:
			# you will need to implement your own FM demodulation based on:
			# https://www.embedded.com/dsp-tricks-frequency-demodulation-algorithms/
			# see more comments on fmSupportLib.py - take particular notice that
			# you MUST have also "custom" state-saving for your own FM demodulator
			fm_demod, state_phase, state_Q = myDemod(i_ds, q_ds, state_phase,state_Q)
		fm_up = [0 for i in range(len(fm_demod)*up)]
		for i in range(len(fm_demod)):
			fm_up[i*up] = fm_demod[i]
            
		# extract the mono audio data through filtering
		if il_vs_th == 0:
			# to be updated by you during the in-lab session based on lfilter
			# same principle as for i_filt or q_filt (but different arguments)
			audio_filt, state_audio_lpf_100k = signal.lfilter(audio_coeff, 1.0, \
						fm_demod,
						zi=state_audio_lpf_100k)
		else:
			# to be updated by you for the takehome exercise
			# with your own code for BLOCK convolution
			audio_filt,state_audio_lpf_100k = block_processing_with_state(audio_coeff,fm_up,state_audio_lpf_100k)


		# downsample audio data
		# to be updated by you during in-lab (same code for takehome)
		audio_block = audio_filt[::audio_decim]
		print(audio_block[0:100])
		# concatenate the most recently processed audio_block
		# to the previous blocks stored already in audio_data
		#
		audio_data = np.concatenate((audio_data, audio_block))
		#

		# to save runtime select the range of blocks to log data
		# this includes both saving binary files as well plotting PSD
		# below we assume we want to plot for graphs for blocks 10 and 11
		if block_count >= 10 and block_count < 12:

			# plot PSD of selected block after FM demodulation
			ax0.clear()
			fmPlotPSD(ax0, fm_demod, (rf_Fs/rf_decim)/1e3, subfig_height[0], \
					'Demodulated FM (block ' + str(block_count) + ')')
			# output binary file name (where samples are written from Python)
			fm_demod_fname = "../data/fm_demod_" + str(block_count) + ".bin"
			# create binary file where each sample is a 32-bit float
			fm_demod.astype('float32').tofile(fm_demod_fname)

			# plot PSD of selected block after extracting mono audio
			# ... change as needed
			ax1.clear()
			fmPlotPSD(ax1, audio_filt, (rf_Fs/rf_decim)/1e3, subfig_height[1], 'Extracted Mono')
			# plot PSD of selected block after downsampling mono audio
			ax2.clear()
			fmPlotPSD(ax2, audio_block, audio_Fs/1e3, subfig_height[2], 'Downsampled Mono Audio')
			# ... change as needed

			# save figure to file
			fig.savefig("../data/myfmMonoBlock" + str(block_count) + ".png")

		block_count += 1

	print('Finished processing all the blocks from the recorded I/Q samples')

	# write audio data to file
	out_fname = "../data/fmMonoBlock.wav"
	wavfile.write(out_fname, int(audio_Fs), np.int16((audio_data/2)*32767))
	print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")

	# uncomment assuming you wish to show some plots
	plt.show()
