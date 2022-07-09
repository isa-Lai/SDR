import numpy as np
import math

def fmPll_RDS(pllIn, freq, Fs, state, ncoScale = 1.0, phaseAdjust = 0.0,normBandwidth = 0.01):
    """
	pllIn 	 		array of floats
					input signal to the PLL (assume known frequency)

	freq 			float
					reference frequency to which the PLL locks

	Fs  			float
					sampling rate for the input/output signals

	ncoScale		float
					frequency scale factor for the NCO output

	phaseAdjust		float
					phase adjust to be added to the NCO output only

	normBandwidth	float
					normalized bandwidth for the loop filter
					(relative to the sampling rate)

	state 			to be added

	"""

	# scale factors for proportional/integrator terms
	# these scale factors were derived assuming the following:
	# damping factor of 0.707 (1 over square root of 2)
	# there is no oscillator gain and no phase detector gain
    Cp = 2.666
    Ci = 3.555
    Kp = (normBandwidth)*Cp
    Ki = (normBandwidth*normBandwidth)*Ci

    ncoOut = np.empty(len(pllIn)+1)
    ncoOutQ = np.empty(len(pllIn)+1)

    integrator = state[0]
    phaseEst = state[1]
    feedbackI = state[2]
    feedbackQ = state[3]
    ncoOut[0] = state[4]
    trigOffset = state[5]

    for k in range(len(pllIn)):
        # phase detector
        errorI = pllIn[k] * (+feedbackI)  # complex conjugate of the
        errorQ = pllIn[k] * (-feedbackQ)  # feedback complex exponential
        # four-quadrant arctangent discriminator for phase error detection
        errorD = math.atan2(errorQ, errorI)
        # loop filter
        integrator = integrator + Ki*errorD
        # update phase estimate
        phaseEst = phaseEst + Kp*errorD + integrator
        # internal oscillator
        trigArg = 2*math.pi*(freq/Fs)*(trigOffset+k+1) + phaseEst
        feedbackI = math.cos(trigArg)
        feedbackQ = math.sin(trigArg)
        ncoOut[k+1] = math.cos(trigArg*ncoScale + phaseAdjust)
        ncoOutQ[k+1] = math.sin(trigArg*ncoScale + phaseAdjust)

    state[0] = integrator
    state[1] = phaseEst
    state[2] = feedbackI
    state[3] = feedbackQ
    state[4] = ncoOut[-1]
    state[5] = trigOffset + len(pllIn)

    return ncoOut, ncoOutQ,state
    # for stereo only the in-phase NCO component should be returned
	# for block processing you should also return the state
    # for RDS add also the quadrature NCO component to the output

if __name__ == "__main__":
    pass
