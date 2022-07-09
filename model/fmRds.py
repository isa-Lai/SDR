
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
from fmPll import fmPll
from fmRRC import impulseResponseRootRaisedCosine
from fmSupportLib import myLowPassFilter,convolveFIR,rational_resampler
from fm_Pll_RDS import fmPll_RDS
# for take-home add your funcions

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
audio_taps = 101#... change as you see fit

up = 17
down = 64
Fs = 240e3
SPS = 27
RRC_Fs = SPS*2375


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


    # FM demodulator (check the library)
    fm_demod, dummy = fmDemodArctan(i_ds, q_ds)
    # we use a dummy because there is no state for this single-pass model


    #RDS channel extraction
    rds_signal_bps = signal.firwin(rf_taps, [54e3/(Fs/2), 60e3/(Fs/2)], window=('hann'), pass_zero="bandpass")
    rds_signal = signal.lfilter(rds_signal_bps,1.0,fm_demod)

    #RDS carrier recovery
    rds_squaring = np.square(rds_signal)
    rds_carrier_bps = signal.firwin(rf_taps,[113.5e3/(Fs/2),114.5e3/(Fs/2)],window = ('hann'),pass_zero = "bandpass")
    rds_carrier = signal.lfilter(rds_carrier_bps,1.0,rds_squaring)

    #pll
    state = [0.0, 0.0, 1.0, 0.0, 1.0, 0.0]
    ncoScale = 1/2
    phase_adjust = math.pi/3.4-math.pi/1.5
    rds_pll_i, rds_pll_q, state = fmPll_RDS(rds_carrier,114e3,Fs,state,ncoScale,phase_adjust,0.001)

    #RDS demodulation
    #mixer
    rds_cha_i = np.zeros(len(rds_signal))
    rds_cha_q = np.zeros(len(rds_signal))
    for i in range(len(rds_signal)):
        rds_cha_i[i] = 2* rds_signal[i]*rds_pll_i[i]
        rds_cha_q[i] = 2*rds_signal[i]*rds_pll_q[i]

    rds_demod_lpf = signal.firwin(rf_taps,3e3/(Fs),window=('hann'))
    rds_lpf_i = signal.lfilter(rds_demod_lpf,1.0,rds_cha_i)
    rds_lpf_q = signal.lfilter(rds_demod_lpf,1.0,rds_cha_q)

    #rational resampler
    anti_img_coeff = signal.firwin(rf_taps,(RRC_Fs/2)/((Fs*up)/2),window=('hann'))
    resample_i, resample_q = rational_resampler(rds_lpf_i,rds_lpf_q,up, down,anti_img_coeff)

    #root_raised cosine filter
    rds_rrc_filter = impulseResponseRootRaisedCosine(RRC_Fs,151)
    rds_rrc_i = signal.lfilter(rds_rrc_filter,1.0,resample_i)
    rds_rrc_q = signal.lfilter(rds_rrc_filter,1.0,resample_q)

    #clock and data recovery
    offset = 0
    max = np.max(rds_rrc_i[0:SPS])
    for i in range(SPS):
        if(rds_rrc_i[i] == max):
            offset = i
            break


    demod_i = rds_rrc_i[offset::SPS]
    demod_q = rds_rrc_q[offset::SPS]

     #Plot scatter plots in order to tune the PLL
    fig, (p_adjust1) = plt.subplots(nrows=1)
    fig.subplots_adjust(hspace = 1.0)
    p_adjust1.scatter(demod_i, demod_q, s=10)
    p_adjust1.set_ylim(-1.5, 1.5)
    plt.show()

    #plt.plot(10*pre_Pll_rds[10180:10200], c='b')
    #plt.plot(post_Pll[10180:10200], c = 'r')
    #plt.show()
    plt.plot(rds_rrc_i[10000:10240], c = 'r')
    plt.plot(rds_rrc_q[10000:10240], c = 'b')
    plt.show()



    # during FM transmission audio samples in the mono channel will contain
    # the sum of the left and right audio channels; hence, we first
    # divide by two the audio sample value and then we rescale to fit
    # in the range offered by 16-bit signed int representation
    print("Written audio samples to \"" + out_fname + "\" in signed 16-bit format")
