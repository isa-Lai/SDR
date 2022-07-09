import sys, math
import numpy as np
from scipy import signal

sample_rate_table = [2400, 2880, 2304, 1920, 1440, 1152, 960]

if __name__ == "__main__":

	if len(sys.argv[0:]) == 1:
		print('\nRun it as follows:\n')
		print('python fmRateChange.py <inputFile> <outFsID> <inFsID>')
		print('\t<inputFile> holds raw I/Q samples (interleaved/8-bit unsigned)')
		print('\t<outFsID> ID for the output (target) sample rate (defaults to 0)')
		print('\t<inFsID>  ID for the input (source) sample rate (defaults to 1)\n')
		print('Sample rate IDs are:')
		print('\t 0 - 2.4   Msamples/sec (defaut input sample rate)')
		print('\t 1 - 2.88  Msamples/sec (defaut output sample rate)')
		print('\t 2 - 2.304 Msamples/sec')
		print('\t 3 - 1.92  Msamples/sec')
		print('\t 4 - 1.44  Msamples/sec')
		print('\t 5 - 1.152 Msamples/sec')
		print('\t 6 - 0.96  Msamples/sec\n')
		exit(1)
	elif len(sys.argv[1:]) >= 1:
		inFsID = 0
		outFsID = 1
		in_fname = sys.argv[1]
		if len(sys.argv[2:]) >= 1:
			outFsID = int(sys.argv[2])
			if len(sys.argv[3:]) == 1:
				inFsID = int(sys.argv[3])
		out_fname = in_fname.partition('.')[0] + '_' + str(sample_rate_table[outFsID]) + '.raw'

	Fs_in = sample_rate_table[inFsID] * 1e3
	Fs_out = sample_rate_table[outFsID]  * 1e3

	# print(out_fname, Fs_out, Fs_in)

	raw_data = np.fromfile(in_fname, dtype='uint8')
	# IQ data is normalized between -1 and +1
	iq_data = (raw_data - 128.0)/128.0
	print("Read raw RF data from \"" + in_fname + "\" in unsigned 8-bit format")

	expand = int(Fs_out) / np.gcd(int(Fs_in), int(Fs_out))
	decim = int(Fs_in) / np.gcd(int(Fs_in), int(Fs_out))

	resampled_i = signal.resample_poly(iq_data[0::2], expand, decim)
	resampled_q = signal.resample_poly(iq_data[1::2], expand, decim)

	# reformat resampled IQ data as 8-bit unsigned
	out_data = np.empty(2*len(resampled_i), dtype='uint8')
	for k in range(len(resampled_i)):
		out_data[2*k] = 128+int(resampled_i[k]*127)
		out_data[2*k+1] = 128+int(resampled_q[k]*127)

	# write resampled IQ data as 8-bit unsigned
	out_data.astype('uint8').tofile(out_fname)
	print("Written resampled RF data from \"" + out_fname + "\" in unsigned 8-bit format")
