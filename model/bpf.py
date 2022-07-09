import matplotlib.animation as animation
import matplotlib.pyplot as plt
from scipy.io import wavfile
from scipy import signal
import numpy as np
import math
import sys
from fmSupportLib import fmDemodArctan, fmPlotPSD
def myLowPassFilter(Fb, Fe, Fs, Ntaps):
	h = np.zeros(Ntaps)
	norm_c = ((Fb+Fe)/2.0)/(Fs/2.0)
	norm_f = (Fe-Fb)/(Fs/2.0);
	for i in range(Ntaps):
		if (i == (Ntaps -1)/2):
			h[i] = norm_f
		else:
			h[i] = norm_f*np.sin(np.pi*norm_f/2.0*(i-(Ntaps-1)/2))/(np.pi*norm_f/2.0*(i-(Ntaps-1)/2))
		h[i] = h[i] *(np.sin(i*np.pi/Ntaps))**2 *(np.cos(i*np.pi*norm_c))
	return h

if __name__ == "__main__":
	h = myLowPassFilter(22e3,54e3,2.4e5,151);
	print(h)
	y = signal.butter(4, [22e3, 54e3], btype='band',fs = 2.4e5)
	
	subfig_height = np.array([0.8, 2, 1.6]) # relative heights of the subfigures
	plt.rc('figure', figsize=(7.5, 7.5))	# the size of the entire figure
	fig, (ax0, ax1, ax2) = plt.subplots(nrows=3, gridspec_kw={'height_ratios': subfig_height})
	fig.subplots_adjust(hspace = .6)
	fmPlotPSD(ax0, h , 2.4e5, subfig_height[1], 'bpf')
	plt.show()