/*
Comp Eng 3DY4 (Computer Systems Integration Project)
Copyright by Nicola Nicolici
Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "dy4.h"
#include "filter.h"
#include "fourier.h"
#include "genfunc.h"
#include "iofunc.h"
#include "logfunc.h"
#include <cmath>
#include <chrono>
#include <algorithm>
#include <atomic>
#include <thread>

//-------------------Global parameter-----------------//
int RFFs, IFFs, AudioFs, rf_decim, audio_decim, if_up;
int taps = 101;
int RDS_if_up , RDS_decim , SPS;
int rf_ele_size,audio_ele_size,block_size;
int queue_element = 8;

void fmDemodArctan(std::vector<float> &fm, float &prev_I, float &prev_Q, const std::vector<float> &i,
	const std::vector<float> &q){
	fm.resize(i.size(),0.0);
	for (unsigned int k = 0; k<i.size();k++){
		float temp = std::pow(i[k],2.0)+std::pow(q[k],2.0);
		if(temp == 0)
			fm[k] = 0.0;
		else{
			fm[k] = (i[k]*(q[k]-prev_Q) - q[k]*(i[k]-prev_I))/temp;
		}
		prev_I = i[k];
		prev_Q = q[k];
	}
}

void RF(int mode, std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset, std::atomic<int>& read_offset_RDS)
{
	std::vector<float> rf_coeff;
	impulseResponseLPF(RFFs, 100e3, taps, rf_coeff);
	for (unsigned int block_id=0; ; block_id++)
	{
		// read raw data
		std::vector<float> block_data(block_size);
		readStdinBlockData(block_size, block_id, block_data);
		if((std::cin.rdstate()) != 0){
			std::cerr << "End of input stream readched" << std::endl;
			exit(1);
		}
	//		std::cerr << "Read block " << block_id <<std::endl;
	//------------------------RF front-end block-------------------
		//extract i and q
		std::vector<float> i_data, q_data, i_flit , q_flit, i_ds , q_ds;
		i_data.reserve(block_data.size()/2);
		q_data.reserve(block_data.size()/2);
		for (unsigned int i = 0 ; i <block_data.size(); i++){
			if(i%2 >0){
				q_data.push_back(block_data[i]);
			}
			else{
				i_data.push_back(block_data[i]);
			}
		}
		static std::vector<float> state_i_lpf_100k;
		static std::vector<float> state_q_lpf_100k;
		state_i_lpf_100k.resize(taps-1,0.0);
		state_q_lpf_100k.resize(taps-1,0.0);

		//convolution with downsampling
		conv_ds(i_flit, i_data,q_flit, q_data, rf_coeff, rf_decim, state_i_lpf_100k, state_q_lpf_100k);
    //conv_ds(i_flit, i_data, rf_coeff, rf_decim, state_i_lpf_100k);
		//conv_ds(q_flit, q_data, rf_coeff, rf_decim, state_q_lpf_100k);

		//demodulation
		static float state_I=0.0, state_Q=0.0;
		std::vector<float> fm_demod;
		fmDemodArctan(fm_demod, state_I, state_Q, i_flit, q_flit);

		//copy to queue
		if(mode == 1 || mode == 3)
			while(write_offset.load() >= (read_offset.load()+queue_element));
		else while(write_offset.load() >= (read_offset.load()+queue_element) ||  write_offset.load() >= (read_offset_RDS.load()+queue_element) );
			//std::this_thread::sleep_for(std::chrono::milliseconds(10));
		std::vector<float>::difference_type address_offset = (write_offset.load() % queue_element)*rf_ele_size;
		std::copy_n(fm_demod.begin(), fm_demod.size(), rf_queue.begin()+address_offset);
		write_offset.fetch_add(1);

	}
}

void MonoStereo(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset,std::vector<float>&audio_queue, std::atomic<int>& audio_write_offset, std::atomic<int>& audio_read_offset)
{
//---------------fix variable----------------
	short int taps_au = (taps-1) * if_up +1;
	std::vector<float> audio_coeff(taps_au);
  impulseResponseLPF_amp(IFFs * if_up, 16e3, taps_au, if_up, audio_coeff);
  std::vector<float> stereo_coeff(taps_au);
  impulseResponseLPF_amp(IFFs * if_up, 15e3, taps_au, if_up, stereo_coeff);
	std::vector<float> stereo_carrier_BPF;
	impulseResponseBPF(IFFs,18.5e3 ,19.5e3, taps, stereo_carrier_BPF);
	std::vector<float> stereo_output_BPF;
	impulseResponseBPF(IFFs,22e3 ,54e3, taps, stereo_output_BPF);
	//std::vector<float> nco;
	std::vector<float> fm_demod(rf_ele_size);
  std::vector<float> mono_allpass;
  static std::vector<float> state_mono_allpass((taps-1)/2);

  std::vector<float> mono;
  static std::vector<float> state_audio(taps_au-1);

  std::vector<float> stereo_carrier;
  static std::vector<float> sc_state_audio(taps-1);

  static float integrator = 0.0 , phaseEst = 0.0, rigOffset = 0.0 ,nextncoOut = 1.0; //state input
	static float feedbackI = 1.0 , feedbackQ = 0.0;
	std::vector<float> ncoOut;

  std::vector<float> stereo_output;
	static std::vector<float> so_state_audio(taps-1);

  std::vector<float> stereo_cha(fm_demod.size());
  std::vector<float> stereo(fm_demod.size()*if_up/audio_decim);
  static std::vector<float> stereo_state_audio(taps_au-1);

  std::vector<float> audio(stereo.size() * 2);

  std::vector<short int> audio_data(audio.size());

	unsigned int fm_size = fm_demod.size();
	unsigned int audio_halfsize = audio.size()/2;
	while(1)
	{
//--------------------------Copy Queue-------------------------------
		while(write_offset.load() <= read_offset.load());
			//std::this_thread::sleep_for(std::chrono::milliseconds(10));
		std::vector<float>::difference_type address_offset = (read_offset.load() % queue_element)*rf_ele_size;
		std::copy_n(rf_queue.begin()+address_offset, fm_demod.size(),fm_demod.begin());
		read_offset.fetch_add(1);
//--------------------------Mono-------------------------------
		allpass(fm_demod, state_mono_allpass, mono_allpass);

		// conv_ds(audio, fm_demod, audio_coeff,audio_decim, state_audio_lpf);

		resampling(mono, mono_allpass, audio_coeff, audio_decim, if_up, state_audio);


//--------------------------Stereo-------------------------------

//carrier
	convolution_with_state(stereo_carrier, sc_state_audio, fm_demod,stereo_carrier_BPF );

	//pll
	fmPll(feedbackI,feedbackQ,integrator, phaseEst, rigOffset, nextncoOut, ncoOut, stereo_carrier, 19e3,IFFs,2 );//   phaseAdjust,  normBandwidth not added yet

//output
	convolution_with_state(stereo_output, so_state_audio, fm_demod,stereo_output_BPF );

//mixer
	for (unsigned int i = 0; i<fm_size;i++)
	{
		stereo_cha[i] = 2* stereo_output[i] * ncoOut[i];
	}

//resample
	resampling(stereo, stereo_cha, stereo_coeff, audio_decim, if_up, stereo_state_audio); //reuse audio_coeff
//combine
	//left audio
	//std::vector<float> audio_l(stereo.size());
	//std::vector<float> audio_r(stereo.size());
  /*
	for (unsigned int i = 0; i<audio_l.size();i++)
	{
		audio_l[i] = (mono[i] +stereo[i])/2;
		audio_r[i] = (mono[i] -stereo[i])/2;
	}*/

	for (unsigned int i = 0; i < audio_halfsize; i++)
	{
			audio[2*i] = (mono[i] +stereo[i])/2;
			audio[2*i+1] = (mono[i] -stereo[i])/2;

}
/*
	nco.insert(nco.end(),ncoOut.begin(), ncoOut.end());
	std::vector<float> freq, psd;
	std::vector<float> freq1, psd1;
	std::vector<float> vector_index;
	genIndexVector(vector_index, nco.size());
	logVector("demod_time", vector_index, nco);
	estimatePSD(fm_demod,IFFs,freq,psd ); //2.4e3 kHz / 10 downsample
	logVector("fmdemod", freq, psd);
	estimatePSD(stereo_output,IFFs,freq,psd ); //2.4e3 kHz / 10 downsample
	logVector("mono", freq, psd);
	estimatePSD(stereo_cha,IFFs,freq,psd ); //2.4e3 kHz / 10 downsample
	logVector("stereo_cha", freq, psd);
	estimatePSD(stereo,AudioFs,freq1,psd1 ); //2.4e3 kHz / 10 downsample
	logVector("stereo", freq1, psd1);
	estimatePSD(audio_l,AudioFs,freq1,psd1 ); //2.4e3 kHz / 10 downsample
	logVector("audio_l", freq1,psd1 );
*/
//--------------------------Write to Aplay-------------------------------

  for (uint k=0; k<audio.size();k++){
    if(std::isnan(audio[k])) audio_data[k] = 0;
    else audio_data[k] = static_cast<short int>(audio[k] * 16384);
  }
  fwrite(&audio_data[0], sizeof(short int), audio_data.size(),stdout);
/*
		//copy to queue
		while(audio_write_offset.load() >= (audio_read_offset.load()+queue_element));
			//std::this_thread::sleep_for(std::chrono::milliseconds(10));
		std::vector<float>::difference_type audio_address_offset = (audio_write_offset.load() % queue_element)*audio_ele_size;
		std::copy_n(audio.begin(), audio.size(), audio_queue.begin()+audio_address_offset);
		audio_write_offset.fetch_add(1);
*/
	}
}


void RDS(std::vector<float>&rf_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset)
{

	//---------------------Fixed variable-------------------
	float RDS_Fs  = SPS * 2375;
	std::vector<float> fm_demod(rf_ele_size);
	unsigned int if_size = fm_demod.size();

	std::vector<float> RDS_output60_BPF;
	impulseResponseBPF(IFFs,54e3 ,60e3, taps, RDS_output60_BPF);
	std::vector<float> RDS_output114_BPF;
	impulseResponseBPF(IFFs,113.5e3 ,114.5e3, taps, RDS_output114_BPF);

	short int taps_au = (taps-1) * RDS_if_up +1;
  std::vector<float> RDS_lp_coeff(taps_au);
  impulseResponseLPF_amp(IFFs * RDS_if_up, 3e3, taps_au, RDS_if_up, RDS_lp_coeff);
	std::vector<float> RDS_RRC;
	fmRRC(RDS_Fs, taps,RDS_RRC );

	std::vector<float> RDS_square_output(rf_ele_size);

	std::vector<float> RDS_output;
	static std::vector<float> RDS_state_audio(taps-1);

	std::vector<float> RDS_carrier_output;
	static std::vector<float> RDS_carrier_state_output(taps-1);

	std::vector<float> RDS_allpass;
	static std::vector<float> state_RDS_allpass((taps-1)/2);

	std::vector<float> RDS_resample_I;
	static std::vector<float> state_RDS_lowpass_I(taps_au-1);
	std::vector<float> RDS_resample_Q;
	static std::vector<float> state_RDS_lowpass_Q(taps_au-1);

	std::vector<float> RDS_RRC_I;
	static std::vector<float> state_RDS_RRC_I(taps-1);
	std::vector<float> RDS_RRC_Q;
	static std::vector<float> state_RDS_RRC_Q(taps-1);

	std::vector<float> RDS_cha_I(fm_demod.size());
	std::vector<float> RDS_cha_Q(fm_demod.size());

	static float integrator = 0.0 , phaseEst = 0.0, trigOffset = 0.0 , phaseAdjust = -1.0*PI, normBandwidth = 0.0; //state input
	std::vector<float> ncoOut_I, ncoOut_Q;

	float  feedbackI = 1.0, feedbackQ = 0.0, nextIOut = 1.0, nextQOut = 0 ;

//---------------------For data recovery-------------------
	int isFirstBlock = 1;
	int sample_offset = 0;
	int isStartOdd = 0; //0 for start even. 1 for start odd

	std::vector<float> symbols_I(rf_ele_size*RDS_if_up/RDS_decim/SPS);
	std::vector<float> symbols_Q(rf_ele_size*RDS_if_up/RDS_decim/SPS);

	std::vector<int> bits(symbols_I.size()/2+1);
	float left_symbol = 1;
	int bits_size;

	std::vector<int> bits_diff(symbols_I.size()/2+1);
	int left_bit = 1;

	std::vector<std::vector<int>> H
	{
		{1,0,0,0,0,0,0,0,0,0},
		{0,1,0,0,0,0,0,0,0,0},
		{0,0,1,0,0,0,0,0,0,0},
		{0,0,0,1,0,0,0,0,0,0},
		{0,0,0,0,1,0,0,0,0,0},
		{0,0,0,0,0,1,0,0,0,0},
		{0,0,0,0,0,0,1,0,0,0},
		{0,0,0,0,0,0,0,1,0,0},
		{0,0,0,0,0,0,0,0,1,0},
		{0,0,0,0,0,0,0,0,0,1},
		{1,0,1,1,0,1,1,1,0,0},
		{0,1,0,1,1,0,1,1,1,0},
		{0,0,1,0,1,1,0,1,1,1},
		{1,0,1,0,0,0,0,1,1,1},
		{1,1,1,0,0,1,1,1,1,1},
		{1,1,0,0,0,1,0,0,1,1},
		{1,1,0,1,0,1,0,1,0,1},
		{1,1,0,1,1,1,0,1,1,0},
		{0,1,1,0,1,1,1,0,1,1},
		{1,0,0,0,0,0,0,0,0,1},
		{1,1,1,1,0,1,1,1,0,0},
		{0,1,1,1,1,0,1,1,1,0},
		{0,0,1,1,1,1,0,1,1,1},
		{1,0,1,0,1,0,0,1,1,1},
		{1,1,1,0,0,0,1,1,1,1},
		{1,1,0,0,0,1,1,0,1,1}
	};


	std::vector<int> syn_A{1,1,1,1,0,1,1,0,0,0};
	std::vector<int> syn_B{1,1,1,1,0,1,0,1,0,0};
	std::vector<int> syn_C{1,0,0,1,0,1,1,1,0,0};
	std::vector<int> syn_Cp{1,1,1,1,0,0,1,1,0,0};
	std::vector<int> syn_D{1,0,0,1,0,1,1,0,0,0};
	std::vector<char> hex{'0','1','2','3','4','5','6','7','8','9','a','b','c','d','e','f'};

	int prevbitLen = 0;
	int blockPosition = 0;
	std::vector<int> prevBits;
	std::vector<int> result_sync(10);
	std::vector<int> block;

//---------------------demode message-------------------
		int nextBlock = 0;
		int haveSync = 0;
		int versionCode;
		int decodeControl = 0;
		std::string PS = "        ";

		const char *PTY[32] = {"None","News","Affairs","Info","Sport","Educated","Drama","Culture","Science","Varied","Pop M","Rock M","Easy M",
					"Light M","Classics","Other M","Weather","Finance","Children","Social","Religion","Phone In","Travel","Leisure","Jazz","Conntry","Nation M","Oldies",
					"Folk M","Document","TEST","Alarm!"};

    //Program Service name
    const std::vector<std::vector<char>> PSN{{' ','0',' ','P'},{' ','1','A','Q'},{' ','2','B','R'},{' ','3','C','S'},{' ','4','D','T'},
          {' ','5','E','U'},{' ','6','F','V'},{'"','7','G','W'},{' ','8','H','X'},{' ','9','I','Y'},{' ',' ','J','Z'},
        {' ',' ','K',' '},{',',' ','L',' '},{'-',' ','M',' '},{'.',' ','N',' '},{'/',' ','O',' '}};

		int syncCout = 0;

	while(1)
	{
//--------------------------Copy Queue-------------------------------
		while(write_offset.load() <= read_offset.load());
			//std::this_thread::sleep_for(std::chrono::milliseconds(10));
		std::vector<float>::difference_type address_offset = (read_offset.load() % queue_element)*rf_ele_size;
		std::copy_n(rf_queue.begin()+address_offset, fm_demod.size(),fm_demod.begin());
		read_offset.fetch_add(1);


	//---------------------RDS Channel Extraction-------------------
	convolution_with_state(RDS_output, RDS_state_audio, fm_demod,RDS_output60_BPF );


	//allpass(RDS_output, state_RDS_allpass, RDS_allpass);
	//---------------------RDS Carrier Recovery-------------------

    //---------Squaring Nonlinearity-----------------
	for(unsigned int i = 0; i< if_size; i++){
    RDS_square_output[i] = RDS_output[i] *RDS_output[i];
	}
	//Carrier
	convolution_with_state(RDS_carrier_output, RDS_carrier_state_output, RDS_square_output,RDS_output114_BPF );
	//PLL&NCO
	RDS_fmPll(feedbackI, feedbackQ, integrator, phaseEst, trigOffset, nextIOut, nextQOut, ncoOut_I, ncoOut_Q, RDS_carrier_output, 114e3,IFFs, 0.5, phaseAdjust,normBandwidth = 0.001);
	//fmPll(feedbackI, feedbackQ, integrator, phaseEst, trigOffset, nextIOut,  ncoOut_I,  RDS_carrier_output, 114e3,IFFs, 0.5, phaseAdjust,normBandwidth = 0.001);
	//Mixer
	for (unsigned int i = 0; i<if_size;i++)
	{
		RDS_cha_I[i] = 2* RDS_output[i] * ncoOut_I[i];
		//RDS_cha_Q[i] = 2* RDS_output[i] * ncoOut_Q[i];
	}
//---------------------RDS Demodulation-------------------
		//resample
	resampling(RDS_resample_I, RDS_cha_I, RDS_lp_coeff, RDS_decim, RDS_if_up, state_RDS_lowpass_I);
	//resampling(RDS_resample_Q, RDS_cha_Q, RDS_lp_coeff, RDS_decim, RDS_if_up, state_RDS_lowpass_Q);

	convolution_with_state(RDS_RRC_I, state_RDS_RRC_I, RDS_resample_I,RDS_RRC );
	//convolution_with_state(RDS_RRC_Q, state_RDS_RRC_Q, RDS_resample_Q,RDS_RRC );

	//clock and data Recovery
	if( isFirstBlock )
	{
		int max = -1;
		for (int i = 0; i<SPS;i++) //to find neex to ignore some first block
		{
			if (std::abs(RDS_RRC_I[i+2*SPS])>max)
			{
				sample_offset = i;
				max = std::abs(RDS_RRC_I[i+2*SPS]);
			}
		}
		//sample_offset = sample_offset%SPS;
		std::cerr <<"Offset: " <<sample_offset << '\n';
	}
	//std::cerr << sample_offset << '\n';
	int index = 0;
	int i = sample_offset;
	for ( ; i<(int)RDS_RRC_I.size();i+=SPS)
	{
		symbols_I[index] = RDS_RRC_I[i];
		//symbols_Q[index] = RDS_RRC_Q[i];
		index++;
	}
	/*
	std::vector<float> vector_index;
	genIndexVector(vector_index, RDS_RRC_I.size());
	logVector("fm_demod", vector_index, symbols_I);
	logVector("fm_demod_q", vector_index, RDS_RRC_I);
	logVector("demod_time", symbols_I, symbols_Q);
	if(nextBlock==3) break;
	nextBlock++;*/
	//next block offset
	//since size is multiple of SPS, change shift is not needed
	//sample_offset = i-RDS_RRC_I.size();

//---------------------RDS Data Process-------------------
	//check if the start point is odd or even
	if(isFirstBlock){
			int num0 = 0, num1 = 0;
			for ( i = 0; i<(int)(symbols_I.size()-2)/2;i++)
			{
				num1 += symbols_I[2*i]*symbols_I[2*i+1] >0 ?1:0;
				num0 += symbols_I[2*i+1]*symbols_I[2*i+2] >0 ?1:0;
			}
			isStartOdd = num1>num0 ? 1:0;
			isFirstBlock = 0;
	}
//convert to bits
	i = 0;
	//deal with left over symbol
	if(isStartOdd)
	{
		bits[0] = left_symbol>symbols_I[0] ? 1:0;
		i = 1;
	}
	//deal with different space
	unsigned int symbol_i_halfsize = (symbols_I.size()+(unsigned int)isStartOdd)/2;
	//if(symbol_i_size != bits.size()){
	bits_size = symbol_i_halfsize;
	for (i = 0; i<bits_size;i++)
	{
		bits[i] = symbols_I[2*i-isStartOdd] > symbols_I[2*i+1-isStartOdd] ? 1:0;
	}
	//check if there is left over
	if((symbols_I.size()+isStartOdd)%2)
	{
		isStartOdd = 1;
		left_symbol = symbols_I[symbols_I.size()-1];
	}
	else isStartOdd = 0;

	//differential coding
	for (i = 0; i<bits_size;i++)
	{
		if(i == 0) bits_diff[i] = bits[i] != left_bit ? 1:0;
		else
		{
				bits_diff[i] = bits[i] != bits[i-1] ? 1:0;
		}
	}
	left_bit = bits[bits_size-1];

//extract message
	//deal with left over bits
	if(prevbitLen!=0)
	{
		blockPosition = 26-prevbitLen;
		std::vector<int> temp(prevBits.begin(), prevBits.end());
		block = temp;
		block.insert(block.end(), bits_diff.begin(), bits_diff.begin()+blockPosition);
	}
	else{
		std::vector<int> temp(bits_diff.begin(), bits_diff.begin()+26);
		block = temp;
		blockPosition = 26;
	}
	std::cerr << "-----------Block Size "<<bits_size << '\n';
	do
	{
		//even already sync, we still check sync for every 30 block;
		if(syncCout>=20){ haveSync = 0; syncCout = 0;}
		//std::cerr << blockPosition << '\n';
		//reset
		for (int k = 0; k<10;k++)
		{
			result_sync[k] = 0;
		}
		if(!haveSync)
		{
			for (int k = 0; k<10;k++)
			{
				for (int l = 0; l<26;l++)
				{
					//tempProd = block[l]*H[l][k];
					result_sync[k] = result_sync[k]^(block[l]*H[l][k]);
					//result_sync[k] = syn_A[k]; //for debug
				}
			}
		}
		//for (int i = 0; i<26;i++) std::cerr << block[i];
		//std::cerr << "" << '\n';
		//check A
		if((haveSync&&nextBlock==0) || result_sync == syn_A)
		{
			std::cerr << "\nBlock A at" << blockPosition -26<< '\n';
			char PICode[4];
			int tempDec;
			for (int m = 0; m<16; m+=4){
					tempDec = block[m]*8 + block[m+1]*4 + block[m+2]*2+ block[m+3]*1;
					PICode[m/4] = hex[tempDec];
			}
			std::cerr << "PI: " << PICode[0] << PICode[1] << PICode[2] << PICode[3] <<   '\n';
			haveSync = 1;
			nextBlock = 1;
		}
		//check b
		else if((haveSync&&nextBlock==1) || result_sync == syn_B)
		{
			std::cerr << "\nBlock B at" << blockPosition -26<< '\n';
			int GroupNum = block[0]*8 + block[1]*4 + block[2]*2+ block[3]*1;
			char GroupChar = block[4]?'B':'A';
			nextBlock =  block[4]?3:2;
			std::cerr << "GroupType "<<GroupNum<<GroupChar << '\n';
			int pty = block[6]*16 +block[7]*8 + block[8]*4 + block[9]*2+ block[10]*1;
			std::cerr << "Program Type: "<<PTY[pty] << '\n';
			decodeControl = block[14]*2+ block[15]*1;
			std::cerr << "Decode Control: "<< decodeControl<< '\n';
			haveSync = 1;

		}
		//check c
		else if((haveSync&&nextBlock==2) || result_sync == syn_C)
		{
			std::cerr << "\nBlock C at" << blockPosition -26<< '\n';
			int altFreq = 0;
			for(int m = 0; m<8;m++){
					altFreq += block[m]*std::pow(2,7-m);
			}
			if(!altFreq) std::cerr << "Alternate Freq No to be Used" << '\n';
			else{
				std::cerr << "Alt Freq(Mhz) "<<87.6+(float)altFreq/10 << '\n';
			}
			haveSync = 1;
			nextBlock = 4;
		}
		//check cp
		else if((haveSync&&nextBlock==3) || result_sync == syn_Cp)
		{
			std::cerr << "\nBlock C' at" << blockPosition -26<< '\n';
			char PICode[4];
			int tempDec;
			for (int m = 0; m<16; m+=4){
					tempDec = block[m]*8 + block[m+1]*4 + block[m+2]*2+ block[m+3]*1;
					PICode[m/4] = hex[tempDec];
			}
			std::cerr << "PI: " << PICode[0] << PICode[1] << PICode[2] << PICode[3] <<   '\n';
			haveSync = 1;
			nextBlock = 4;
		}
		//check d
		else if((haveSync&&nextBlock==4) || result_sync == syn_D)
		{
			std::cerr << "\nBlock D at" << blockPosition -26<< '\n';
			int PSCode[4];
			for (int m = 0; m<16;m+=4)
			{
				PSCode[m/4]  = block[m]*8 + block[m+1]*4 + block[m+2]*2+ block[m+3]*1;
				if(m/4%2==0)
				{
					if(PSCode[m/4]>5) PSCode[m/4] = 0;// for b7-b4 >5, not recoded in table
					else PSCode[m/4]-=2;
				}
			}
			PS[decodeControl*2] = PSN[PSCode[1]][PSCode[0]];
			PS[decodeControl*2+1] = PSN[PSCode[3]][PSCode[2]];
			std::cerr << "Program Service: "<<PS<< '\n';
			haveSync = 1;
			nextBlock = 0;
		}
		//move frame
		if(haveSync)
			{
				//move to next block
				std::vector<int> temp(bits_diff.begin()+blockPosition, bits_diff.begin()+blockPosition+26);
				blockPosition+=26;
				block = temp;
				//haveSync = 0;
				syncCout++;
			}
			else{
			block.erase(block.begin());
			block.push_back(bits_diff[blockPosition]);
			blockPosition++;
			}

	}
	while(blockPosition<bits_size);
	//deal with left over bits
	std::vector<int> temp(bits_diff.begin()+blockPosition-26, bits_diff.begin()+bits_size);
	prevBits = temp;
	prevbitLen = prevBits.size();

	}

}

void AudioWrite(std::vector<float>&audio_queue, std::atomic<int>& write_offset, std::atomic<int>& read_offset)
{
	std::vector<float> audio(audio_ele_size);
	while(1)
	{
		//copy queue
		while(write_offset.load() <= read_offset.load());
		//std::this_thread::sleep_for(std::chrono::milliseconds(10));
		std::vector<float>::difference_type address_offset = (read_offset.load() % queue_element)*audio_ele_size;
		std::copy_n(audio_queue.begin()+address_offset, audio.size(),audio.begin());
		read_offset.fetch_add(1);
		//------------Send Audio---------------
		std::vector<short int> audio_data(audio.size());
		for (uint k=0; k<audio.size();k++){
			if(std::isnan(audio[k])) audio_data[k] = 0;
			else audio_data[k] = static_cast<short int>(audio[k] * 16384);
		}
		fwrite(&audio_data[0], sizeof(short int), audio_data.size(),stdout);

		//ending thread if necessary
	}
}

int main(int argc,char * argv[])
{
//-------------------------read arguement----------------------
	//assume the first input will be the mode
	int mode = 0;

	if(argc<2)
	{
		std::cerr << "Operating in default mode 0" << std::endl;
	} else if (argc == 2) {
		mode = atoi(argv[1]);
		if (mode > 3) {
			std::cerr << "Wrong mode" << mode << std::endl;
			exit(1);
		}
	} else {
		std::cerr <<"Usage: "<< argv[0] << std::endl;
		std::cerr << "or" << std::endl;
		std::cerr << "Usage: " << argv[0] << " <mode>" << std::endl;
		std::cerr << "\t\t <mode> is a value from 0 to 3" << std::endl;
		exit(1);
	}
//	int mode = (int)(*argv[1]) - 48;
	std::cerr << "Operating in mode "<<mode << '\n';

	//parameter about mode
		switch (mode){
			case 1:
				RFFs = 1.44e6;
				IFFs = 2.88e5;
				AudioFs = 4.8e4;
				rf_decim = 5;
				audio_decim = 6;
				if_up = 1;
				block_size = 200 * rf_decim * audio_decim * 2;
				break;
			case 2:
				RFFs = 2.4e6;
				IFFs = 2.4e5;
				AudioFs = 4.41e4;
				rf_decim = 10;
				RDS_decim = 48;
				audio_decim = 800;
				if_up = 147;
				RDS_if_up = 19;
				//block_size = 50 * rf_decim * audio_decim * 2;
				SPS = 40;
				block_size = 38*2*SPS * rf_decim * RDS_decim * 2/RDS_if_up;
				break;
			case 3:
				RFFs = 1.92e6;
				IFFs = 3.84e5;
				AudioFs = 4.41e4;
				rf_decim = 5;
				audio_decim = 1280;
				if_up = 147;
				block_size = 1 * rf_decim * audio_decim * 2;
				break;
			default:
				//mode 0
				RFFs = 2.4e6;
				IFFs = 2.4e5;
				AudioFs = 4.8e4;
				rf_decim = 10;
				RDS_decim = 640;
				audio_decim = 5;
				if_up = 1;
				RDS_if_up = 171;
				SPS = 27;
				//block_size = 200 * rf_decim * audio_decim * 2;
				block_size = 38*2*SPS * rf_decim * RDS_decim * 2/RDS_if_up;//dont touch this
		}
//---------------------- Raad raw & RF Front End-------------------
	//block_size = 200 * rf_decim * audio_decim * 2;
	rf_ele_size = block_size/2/rf_decim;
	std::vector<float> rf_queue(rf_ele_size*queue_element);
	std::atomic<int> rf_write_offset(0);
	std::atomic<int> path_read_offset(0);
	std::atomic<int> RDS_read_offset(0);

	std::thread rffrontend = std::thread(RF, std::ref(mode), std::ref(rf_queue),std::ref(rf_write_offset),std::ref(path_read_offset),std::ref(RDS_read_offset));
//-----------------------Mono Path---------------------------
	audio_ele_size = block_size/rf_decim/audio_decim*if_up;
	std::vector<float> audio_queue(audio_ele_size*queue_element);
	std::atomic<int> path_write_offset(0);
	std::atomic<int> audio_read_offset(0);
	std::thread audiopath = std::thread(MonoStereo, std::ref(rf_queue),std::ref(rf_write_offset),std::ref(path_read_offset),std::ref(audio_queue),std::ref(path_write_offset),std::ref(audio_read_offset));


	//std::vector<float> down_sampled;
	 //estimatePSD(stereo_carrier, (AudioFs/1e3), 1, down_sampled);
	 //logVector("stereo_carrier", freq, down_sampled);
	//std::cout << "Run: gnuplot -e 'set terminal png size 1024,768' ../data/example.gnuplot > ../data/example.png\n";


//-----------------------RDS Path-------------------------------
	if(mode ==0 || mode == 2)
	{
		std::thread RDSpath = std::thread(RDS, std::ref(rf_queue),std::ref(rf_write_offset),std::ref(RDS_read_offset));
		RDSpath.join();
	}
  else{
    std::cerr << "Mode "<<mode << "does not support RDS\n";
  }

//-------------------Finishing write-----------
	//std::thread writeaudio = std::thread(AudioWrite, std::ref(audio_queue),std::ref(path_write_offset),std::ref(audio_read_offset));

	rffrontend.join();
	audiopath.join();
	//writeaudio.join();
	return 0;
}
