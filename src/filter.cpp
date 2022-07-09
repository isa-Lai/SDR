/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"

void fmRRC(float Fs, unsigned short int num_taps, std::vector<float> &impulseResponseRRC){
	/*
	Root raised cosine (RRC) filtered

	Fs			sampling rate at the output of the
	resampler in the RDS path
					sampling rate must be an integer multiplier
					of 2375
					this integer multiple is the number of
					samples per symbol
	num_taps number of filter taps
	*/
	// duration for each symbol - should NOT be changed for RDS!
	float T_symbol = 1.0/2375.0;

	// roll-off factor(must be greater than 0 and smaller than1)
	// for RDS a vaule in the range of 0.9 is a good trade-off between
	// the excess bandwidth and the size/ duration of ripples in the time-domain
	float beta = 0.90;

	// the RRC inpulse response that will be coumputed in this function
	impulseResponseRRC.clear();
	impulseResponseRRC.resize(num_taps,0.0);

	for(unsigned int k = 0; k<num_taps;k++){
		float t = ((float)k-(float)num_taps/2.0)/Fs;
		if( t == 0.0){
			impulseResponseRRC[k] = 1.0 + beta *(4.0/PI - 1.0);
		}
		else if (t == (0.0-T_symbol)/(4.0*beta) || t== T_symbol/(4.0*beta)){
			impulseResponseRRC[k] = (beta/1.414213562)*(((1.0+2.0/PI)* (sin(PI/(4.0*beta)))) + ((1.0-2.0/PI)*(cos(PI/(4.0*beta)))));
		}
		else{
			impulseResponseRRC[k] = (sin(PI*t*(1.0-beta)/T_symbol)+ 4.0*beta*(t/T_symbol)* cos(PI*t*(1.0+beta)/T_symbol))/ (PI*t*(1.0-(4.0*beta*t/T_symbol)*(4.0*beta*t/T_symbol))/T_symbol);
		}
	}
	//the RRC impulse response to be used by convolution.
}

void fmPll(float &feedbackI , float &feedbackQ, float &integrator , float &phaseEst, float &trigOffset,float &nextncoOut, std::vector<float> &ncoOut, std::vector<float> &pllIn,  float freq, float Fs, float ncoScale , float phaseAdjust, float normBandwidth )
{
	//default parameter declared in .h

	//scale factors for proportional/integrator terms
	//these scale factors were derived assuming the following:
	//damping factor of 0.707 (1 over square root of 2)
	//there is no oscillator gain and no phase detector gain
	float Cp = 2.666;
	float Ci = 3.555;

	//gain for the proportional term
	float Kp = (normBandwidth)*Cp;
	//gain for the integrator term
	float Ki = (normBandwidth*normBandwidth)*Ci;

	//output array for the NCO
	ncoOut.clear(); ncoOut.resize(pllIn.size(), 0.0);

	//initialize internal state
	float trigArg;
	ncoOut[0] = nextncoOut;

	//note: state saving will be needed for block processing
	float errorI, errorQ, errorD;
	for(unsigned int k = 0; k<pllIn.size(); k++)
	{

		//phase detector
		errorI = pllIn[k] * (+feedbackI);  //# complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ);  //# feedback complex exponential

		//four-quadrant arctangent discriminator for phase error detection
		errorD = std::atan2(errorQ, errorI);

		//loop filter
		integrator = integrator + Ki*errorD;

		//update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator;

		//internal oscillator
		trigOffset += 1;
		trigArg = 2*PI*(freq/Fs)*(trigOffset) + phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		if(k+1 == pllIn.size()) nextncoOut = std::cos(trigArg*ncoScale + phaseAdjust);
		else ncoOut[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);

	//for stereo only the in-phase NCO component should be returned
	//for block processing you should also return the state
	//for RDS add also the quadrature NCO component to the output
	}
}

void RDS_fmPll(float &feedbackI , float &feedbackQ, float &integrator , float &phaseEst, float &trigOffset,float &nextIOut, float &nextQOut, std::vector<float> &ncoOut_I, std::vector<float> &ncoOut_Q, std::vector<float> &pllIn,  float freq, float Fs, float ncoScale , float phaseAdjust, float normBandwidth)
{
	//default parameter declared in .h

	//scale factors for proportional/integrator terms
	//these scale factors were derived assuming the following:
	//damping factor of 0.707 (1 over square root of 2)
	//there is no oscillator gain and no phase detector gain
	float Cp = 2.666;
	float Ci = 3.555;

	//gain for the proportional term
	float Kp = (normBandwidth)*Cp;
	//gain for the integrator term
	float Ki = (normBandwidth*normBandwidth)*Ci;

	//output array for the NCO
	ncoOut_I.clear(); ncoOut_I.resize(pllIn.size(), 0.0);
	ncoOut_Q.clear(); ncoOut_Q.resize(pllIn.size(), 0.0);

	//initialize internal state
	float  trigArg;
	ncoOut_I[0] = nextIOut;
	ncoOut_Q[0] = nextQOut;

	//note: state saving will be needed for block processing
	float errorI, errorQ, errorD;
	for(unsigned int k = 0; k<pllIn.size(); k++)
	{

		//phase detector
		errorI = pllIn[k] * (+feedbackI);  //# complex conjugate of the
		errorQ = pllIn[k] * (-feedbackQ);  //# feedback complex exponential

		//four-quadrant arctangent discriminator for phase error detection
		errorD = std::atan2(errorQ, errorI);

		//loop filter
		integrator = integrator + Ki*errorD;

		//update phase estimate
		phaseEst = phaseEst + Kp*errorD + integrator;

		//internal oscillator
		trigOffset += 1;
		trigArg = 2*PI*(freq/Fs)*(trigOffset) + phaseEst;
		feedbackI = std::cos(trigArg);
		feedbackQ = std::sin(trigArg);
		if(k+1 == pllIn.size())
		{
			 nextIOut = std::cos(trigArg*ncoScale + phaseAdjust);
			 nextQOut = std::sin(trigArg*ncoScale + phaseAdjust);
		 }
		else {
		    ncoOut_I[k+1] = std::cos(trigArg*ncoScale + phaseAdjust);
		    ncoOut_Q[k+1] = std::sin(trigArg*ncoScale + phaseAdjust);
		}

	//for stereo only the in-phase NCO component should be returned
	//for block processing you should also return the state
	//for RDS add also the quadrature NCO component to the output
	}
}

// function to compute the impulse response "h" based on the sinc function
void impulseResponseBPF(float Fs, float Fb,float Fe, unsigned short int num_taps, std::vector<float> &h)
{
	// bring your own functionality
  // allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	const float normCenter = ((Fb+Fe)/2.0)/(Fs/2.0);
	const float normPass = (Fe-Fb)/(Fs/2.0);
	for (unsigned short int i = 0; i<num_taps;i++){
		if(i == (num_taps-1)/2.0){
			h[i] = normPass;
		}
		else{
			h[i] = normPass*(std::sin(PI*normPass/2.0*(i-(num_taps-1)/2.0)))/(PI*normPass/2.0*(i-(num_taps-1)/2.0));
		}
		h[i] = h[i] * std::cos(i*PI*normCenter)* std::pow(std::sin(i*PI/static_cast<float>(num_taps)),2.0);
	}
}


// function to compute the impulse response "h" based on the sinc function
void impulseResponseLPF(float Fs, float Fc, unsigned short int num_taps, std::vector<float> &h)
{
	// bring your own functionality
  // allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	const float normCut = Fc/(Fs/2.0);
	for (unsigned short int i = 0; i<num_taps;i++){
		if(i == (num_taps-1)/2.0){
			h[i] = normCut;
		}
		else{
			h[i] = normCut*(std::sin(PI*normCut*(i-(num_taps-1)/2.0)))/(PI*normCut*(i-(num_taps-1)/2.0));
		}
		h[i] = h[i] * std::pow(std::sin(i*PI/static_cast<float>(num_taps)),2.0);
	}
}

void impulseResponseLPF_amp(float Fs, float Fc, unsigned short int num_taps, int amp, std::vector<float> &h)
{
	// allocate memory for the impulse response
	h.clear(); h.resize(num_taps, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	const float normCut = Fc/(Fs/2.0);
	for (unsigned short int i = 0; i<num_taps;i++){
		if(i == (num_taps-1)/2.0){
			h[i] = normCut;
		}
		else{
			h[i] = normCut*(std::sin(PI*normCut*(i-(num_taps-1)/2.0)))/(PI*normCut*(i-(num_taps-1)/2.0));
		}
		h[i] = amp * h[i] * std::pow(std::sin(i*PI/static_cast<float>(num_taps)),2.0);
	}
}

void allpass(std::vector<float> &input_block, std::vector<float> &state_block, std::vector<float> &output_block)
{
	// allocate memory for the impulse response
	output_block.clear();
	output_block = state_block;
	output_block.reserve(input_block.size());
	for(unsigned int n = 0; n<(input_block.size()-state_block.size());n++){
		output_block.push_back(input_block[n]);
	}

	std::vector<float> temp(&input_block[input_block.size()-state_block.size()], &input_block[input_block.size()]);
	state_block = temp;
}

// function to compute the filtered output "y" by doing the convolution
// of the input data "x" with the impulse response "h"
void convolveFIR(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h)
{
	// bring your own functionality
  // allocate memory for the output (filtered) data
	y.clear(); y.resize(x.size()+h.size()-1, 0.0);

	// the rest of the code in this function is to be completed by you
	// based on your understanding and the Python code from the first lab
	for (unsigned int n = 0;n<y.size();n++){
		for (unsigned int k=0;k<h.size();k++){
			if(n>=k && n-k<x.size()){
				y[n] += h[k]*x[n-k];
			}
		}
	}
}

void convolution_with_state(std::vector<float> &yb, std::vector<float> &state, const std::vector<float> &xb, const std::vector<float> &h){
	yb.clear();
	yb.resize(xb.size(), 0.0);

	int n,k,diff;
	int ybsize = (int)yb.size();
	int hsize = (int)h.size();
	for (n=0; n<ybsize; n++){
		for (k=0; k<hsize; k++){
			diff = n-k;
			if (diff >= 0){
				yb[n] = yb[n] + h[k] * xb[diff];
			}else{
				yb[n] = yb[n] + h[k] * state[state.size()+diff];
			}
		}
	}
	std::vector<float> temp(&xb[xb.size()-state.size()], &xb[xb.size()]);
	state = temp;
}
//Note!!!!!!!!!!! this function is only for i and q data filtering
void conv_ds(std::vector<float> &yi, const std::vector<float> &xi, std::vector<float> &yq, const std::vector<float> &xq,const std::vector<float> &h, const int ds, std::vector<float> &statei,std::vector<float> &stateq)
{
	yi.clear();
	yi.resize(xi.size()/ds, 0.0);
	yq.clear();
	yq.resize(xq.size()/ds, 0.0);

	for (int n = 0;n<(int)yi.size();n++){
		for (int k = 0;k<(int)h.size();k++){
			int idx = ds*n-k;
			if( idx >= 0){
				yi[n] += h[k] * xi[idx];
				yq[n] += h[k] * xq[idx];
			}else {
				yi[n] += h[k] * statei[statei.size() - abs(idx)];
				yq[n] += h[k] * stateq[stateq.size() - abs(idx)];
			}
		}
	}
	std::vector<float> tempi(&xi[xi.size()-statei.size()], &xi[xi.size()]);
	statei = tempi;
	std::vector<float> tempq(&xq[xq.size()-stateq.size()], &xq[xq.size()]);
	stateq = tempq;
}
/*
//convolution with downsampling
void conv_ds(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, const int ds, std::vector<float> &state)
{
	y.clear();
	y.resize(x.size()/ds, 0.0);

	for (int n = 0;n<(int)y.size();n++){
		for (int k = 0;k<(int)h.size();k++){
			if((ds*n-k) >= 0){
				y[n] += h[k] * x[(ds*n-k)];
			}else {
				y[n] += h[k] * state[state.size() - abs((ds*n-k))];
			}
		}
	}
	std::vector<float> temp(&x[x.size()-state.size()], &x[x.size()]);
	state = temp;
}
*/
void resampling(std::vector<float> &y, const std::vector<float> &x, const std::vector<float> &h, const int ds,
								const int us, std::vector<float> &state)
{
		y.clear(); y.resize(x.size()*us/ds, 0.0);
		int n,k,idx;
		int ybsize = (int)y.size();
		int hsize = (int)h.size();
		for (n = 0;n<ybsize;n++){
			for (k=(n*ds)%us ;k<hsize; k=k+us){
				idx = (n*ds-k)/us;
				if( idx >= 0 ){
					y[n] += h[k]*x[idx];
				}
				else {
					y[n] += h[k] * state[state.size() +idx];
				}
			}
		}
	std::vector<float> temp(&x[x.size()-state.size()], &x[x.size()]);
	state = temp;
}
