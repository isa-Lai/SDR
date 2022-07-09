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

from fmSupportLib import fmDemodArctan, fmPlotPSD
from fmPll import fmPll
from fmRRC import impulseResponseRootRaisedCosine
from fmSupportLib import rational_resampler
from fm_Pll_RDS import fmPll_RDS
# for take-home add your funcions


rf_Fs = 2.4e6
rf_Fc = 100e3
rf_taps = 101
rf_decim = 10
Fs = 240e3
audio_Fs = 240e3
audio_decim = 3
up = 171
down = 640
SPS = 27
RRC_Fs = SPS*2375
# add other settings for audio, like filter taps, ...
audio_taps = 101
audio_Fc = 16e3
block_size = 1024*rf_decim*audio_decim*10

if __name__ == "__main__":

    # read the raw IQ data from the recorded file
    # IQ data is assumed to be in 8-bits unsigned (and interleaved)
    in_fname = "../data/samples8.raw"
    raw_data = np.fromfile(in_fname, dtype='uint8')
    print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")
    # IQ data is normalized between -1 and +1 in 32-bit float format
    iq_data = (np.float32(raw_data) - 128.0)/128.0
    print("Reformatted raw RF data to 32-bit float format (" + str(iq_data.size * iq_data.itemsize) + " bytes)")

    # select a block_size that is a multiple of KB
    # and a multiple of decimation factors
    #block_size = 1024 * rf_decim * audio_decim * 2
    block_count = 0
    iq_data = iq_data[0:6*block_size]

    rf_coeff = signal.firwin(rf_taps, rf_Fc/(rf_Fs/2), window=('hann'))

    # states needed for continuity in block processing
    state_i_lpf_100k = np.zeros(rf_taps-1)
    state_q_lpf_100k = np.zeros(rf_taps-1)
    state_audio_lpf_100k = np.zeros(audio_taps*up-1)
    state_phase = 0

   #RDS coefficients
    #RDS extraction
    rds_signal_bpf = signal.firwin(rf_taps, [54e3/(audio_Fs/2), 60e3/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
    rds_signal_state = np.zeros(rf_taps-1)

    #RDS_carrier recovery
    rds_carrier_bpf = signal.firwin(rf_taps, [113.5e3/(audio_Fs/2), 114.e3/(audio_Fs/2)], window=('hann'), pass_zero="bandpass")
    rds_carrier_state =np.zeros(rf_taps-1)
    #pll
    freq_centered =114000
    phase_adj = -0.25*math.pi
    state_Pll =[0.0, 0.0, 1.0, 0.0, 1.0, 0.0]

    #3k LPF
    rds_carrier_lpf = signal.firwin(rf_taps, 3e3/(audio_Fs/2), window=('hann'))
    rds_carrier_lpf_state_i = np.zeros(rf_taps-1)
    rds_carrier_lpf_state_i_q = np.zeros(rf_taps-1)

    #RDS demodulation
    #rational resampler
    up = 171
    down = 640
    SPS = 27
    rrc_Fs = SPS*2375
    resample_filt_coeff = signal.firwin(rf_taps, (rrc_Fs/2)/((240000*up)/2), window=('hann'))
    resample_state_i = np.zeros(rf_taps-1)
    resample_state_q = np.zeros(rf_taps-1)
    #Values for RRC
    rrc_taps = 101
    rds_rrc_filt = impulseResponseRootRaisedCosine(rrc_Fs, rrc_taps)
    rds_rrc_state_i = np.zeros(rrc_taps -1)
    rds_rrc_state_q = np.zeros(rrc_taps -1)
    #Values for clock recoverey
    offset = 0

    while (block_count+1)*block_size < len(iq_data):
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
        fm_demod, state_phase = fmDemodArctan(i_ds, q_ds, state_phase)

    # RDS extraction
        rds_signal, rds_signal_state = signal.lfilter(rds_signal_bpf,1.0,fm_demod,zi=rds_signal_state)

        # Carrier Recovery
        #Squaring
        rds_squaring = np.square(rds_signal)
        rds_carrier,rds_carrier_state = signal.lfilter(rds_carrier_bpf,1.0,rds_squaring,zi=rds_carrier_state)

        #pll
        rds_pll_i, rds_pll_q, state_Pll =  fmPll_RDS(rds_carrier, freq_centered, 240000, state_Pll, ncoScale = 0.5, phaseAdjust =phase_adj , normBandwidth = 0.001)

        #RDS demodulation
        # Mixer
        rds_mix_i = 2*np.multiply(rds_signal, rds_pll_i[0:len(rds_signal):1])
        rds_mix_q = 2*np.multiply(rds_signal, rds_pll_q[0:len(rds_signal):1])

        #LPF
        rds_carrier_lpf_i,rds_carrier_lpf_state_i = signal.lfilter(rds_carrier_lpf,1.0, rds_mix_i,zi=rds_carrier_lpf_state_i)
        rds_carrier_lpf_q,rds_carrier_lpf_state_q = signal.lfilter(rds_carrier_lpf,1.0, rds_mix_q,zi=rds_carrier_lpf_state_i_q)
        upsample_i = np.zeros(len(rds_carrier_lpf_i)*up)
        upsample_q = np.zeros(len(rds_carrier_lpf_i)*up)

        #Rational Resampler
        #Upsamples
        for i in range (len(rds_carrier_lpf_i)):
            upsample_i[i*up] = rds_carrier_lpf_i[i]
            upsample_q[i*up] = rds_carrier_lpf_q[i]
        #Downsample
        rds_resample_filt_i,resample_state_i= signal.lfilter(resample_filt_coeff,1.0, upsample_i,zi=resample_state_i)
        rds_resample_filt_q,resample_state_q= signal.lfilter(resample_filt_coeff,1.0, upsample_q,zi=resample_state_q)
        resample_i = rds_resample_filt_i[::down]*up
        resample_q = rds_resample_filt_q[::down]*up

        #RRC
        rrc_rds,rds_rrc_state_i = signal.lfilter(rds_rrc_filt, 1.0, resample_i,zi=rds_rrc_state_i)
        rrc_rds_Q, rds_rrc_state_q= signal.lfilter(rds_rrc_filt, 1.0, resample_q,zi=rds_rrc_state_q)

        #Clock and data recovery

        max = np.max(rrc_rds[0:SPS])
        if block_count == 0:
            for i in range(SPS):
                if(rrc_rds[i] == max):
                    offset = i
                    break
        print('offset: ',offset)

        symbols_I = rrc_rds[offset::SPS]
        symbols_Q = rrc_rds_Q[offset::SPS]
        #block processing the offset for the next block
        for i in range(len(rrc_rds)-SPS,len(rrc_rds),1):
            if(rrc_rds[i] == symbols_I[-1]):
                offset = SPS-i
                break

        #Plotting
        if block_count == 0:
            fig, (p_adjust1) = plt.subplots(nrows=1)
            fig, (rrc) = plt.subplots(nrows=1)
            fig.subplots_adjust(hspace = 1.0)
            p_adjust1.scatter(symbols_I, symbols_Q, s=10)
            p_adjust1.set_ylim(-1.25, 1.25)
            rrc.plot(rrc_rds[0:512], c = 'r')
            rrc.plot(rrc_rds_Q[0:512], c = 'b')

            # set up drawing
            fig, (ax0, ax1, ax2,ax4) = plt.subplots(nrows=4)
            fig.subplots_adjust(hspace = 1.0)

            # PSD after FM demodulation
            ax0.psd(fm_demod, NFFT=512, Fs=(rf_Fs/rf_decim)/1e3)
            ax0.set_ylabel('PSD (db/Hz)')
            ax0.set_title('Demodulated FM')

            # save PSD plots
            ax1.psd(rds_signal, NFFT=512, Fs=(240e3/rf_decim)/1e3)
            ax1.set_ylabel('PSD (db/Hz)')
            ax1.set_title('Pre Pll')

            ax2.psd(rrc_rds_Q, NFFT=512, Fs=(240e3/rf_decim)/1e3)
            ax2.set_ylabel('PSD (db/Hz)')
            ax2.set_title('Post Pll')

            ax4.psd(rrc_rds, NFFT=512, Fs=(rrc_Fs/1e3))
            ax4.set_ylabel('PSD (db/Hz)')
            ax4.set_title('Post RRC Filter')

            ## RDS

            plt.show()

        block_count += 1

# uncomment assuming you wish to show some plots
plt.show()
