// Template for C6713 DSK and AIC23 codec
#define CHIP_6713

#include "dsk6713.h"
#include "dsk6713_aic23.h" //

#include "data_trans.h"
//#include "data_rec.h"
#include "data_dechirp.h"

//#define USE_FFT_LIB
#include "FFT_helper.h"

DSK6713_AIC23_CodecHandle hCodec;							// Codec handle
DSK6713_AIC23_Config config = DSK6713_AIC23_DEFAULTCONFIG;  // Codec configuration with default settings

///////////////////////////////////
int interrupt_count = 0;
double MAX_SIGNAL_VALUE = 0;

void initialize();
short convert_double_to_short_for_output(double max_value, double input_value);

#ifndef USE_FFT_LIB
void prepare_fft_data(COMPLEX* data, int N){
	// need to fill
	int i;

	for (i = 0; i < N; i++){
		if (i > NUM_SAMPLES_PER_CHIRP){
						data[i].real = 0;
						data[i].imag = 0;
					}
		else{
		data[i].real = dechirp_real[i];
		data[i].imag = dechirp_img[i];
		}

	}

}

#else
	
void prepare_fft_data_for_lib(double* data, int N)
{

	int i;

		for (i = 0; i < N*2; i+=2){


			if (i > NUM_SAMPLES_PER_CHIRP){
				data[i] = 0;
				data[i+1] = 0;
			}
			else{
				data[i] = dechirp_real[i/2];
				data[i+1] = dechirp_img[i/2];
			}
		}

}

#endif

Uint32 get_timer_overhead(TIMER_Handle hTimer){
// Timer Setup
	TIMER_configArgs(hTimer, 0x000002C0, 0xFFFFFFFF, 0x00000000);
	Uint32 start = TIMER_getCount(hTimer);
		   start = TIMER_getCount(hTimer);
	Uint32 stop  = TIMER_getCount(hTimer);

	Uint32 overhead = stop - start;
	return overhead;
}

void main()
{
	interrupt_count =0;

	TIMER_Handle hTimer = TIMER_open(TIMER_DEVANY, 0);	
	Uint32 timer_overhead = get_timer_overhead(hTimer);
	
	initialize();
	
#ifndef USE_FFT_LIB
	prepare_w_fft(w_fft, PTS);
	prepare_fft_data(data_fft, PTS);
#else
	//use prepare_
	prepare_w_fft_lib (w_fft_lib, PTS);
	prepare_fft_data_for_lib(data_fft_lib, PTS);
	//need to fill
#endif

	
	Uint32 timer_start = TIMER_getCount(hTimer);
#ifndef USE_FFT_LIB
	FFT(data_fft, w_fft, PTS);
#else


	DSPF_dp_cfftr2(PTS, data_fft_lib, w_fft_lib, 1);
	bit_rev_2(data_fft_lib, PTS);




	//FFT(data_fft_lib, w_fft_lib, PTS);

	// need to fill
#endif
	Uint32 timer_end = TIMER_getCount(hTimer);	
	Uint32 diff = (timer_end - timer_start) - timer_overhead;
	
	int count_to_cpu_cycles = 4;
	printf("Number of cycles for FFT calculation: %d \n", diff*count_to_cpu_cycles);

	TIMER_close(hTimer);
	
	MAX_SIGNAL_VALUE = 0; //for c_int11 to recalculate the max value
	while(1);					// main loop - do nothing but wait for interrupts
}

/** interrupt service routine */
interrupt void c_int11()		
{
	int num_gaps_to_insert = 50;
	
	double sig = 0;
	if (interrupt_count < PTS){
		// Note: do not use pow(), it has some issue
#ifndef USE_FFT_LIB
		sig = sqrt(data_fft[interrupt_count].real*data_fft[interrupt_count].real + data_fft[interrupt_count].imag*data_fft[interrupt_count].imag);		
#else
		sig = sqrt(data_fft_lib[2*interrupt_count]*data_fft_lib[2*interrupt_count] + data_fft_lib[2*interrupt_count+1]*data_fft_lib[2*interrupt_count+1]);		
#endif
	}
	else if (interrupt_count < PTS + num_gaps_to_insert){
		sig = -0.5*MAX_SIGNAL_VALUE;//0;		
	}
	
	interrupt_count++;
	if (interrupt_count >= NUM_SAMPLES + num_gaps_to_insert){
		interrupt_count = 0;
	}
	
	if (sig > MAX_SIGNAL_VALUE){
		MAX_SIGNAL_VALUE = sig;
	}		
	
	short sig_out = convert_double_to_short_for_output(MAX_SIGNAL_VALUE, sig);
	
	/* DO NOT MODIFY ANYTHING BELOW */
	sig_out = -1*sig_out; //Board inverts the signal from the DAC.
	
	union {Uint32 stereo; short channel[2];} audio;	
	audio.channel[1] = sig_out;	// insert left channel
	audio.channel[0] = sig_out;	// insert right channel
	MCBSP_write(DSK6713_AIC23_DATAHANDLE, audio.stereo); // write stereo audio	
}

/* DO NOT MODIFY ANYTHING BELOW */
void vectors();         // external function

short convert_double_to_short_for_output(double max_value, double input_value)
{	
	short MAX_POSITIVE_VALUE = 32700;//32767
	
	if (input_value > max_value)
		input_value = max_value;
	
	short result = (short)(MAX_POSITIVE_VALUE/max_value*input_value);
	
	return result;
}

void initialize()
{
	// initialize dsp
	DSK6713_init();		//call BSL to init DSK-EMIF,PLL)

	DSK6713_LED_init();

	// initialize AIC23 peripheral
	hCodec = DSK6713_AIC23_openCodec(0, &config);	//handle(pointer) to codec
	DSK6713_AIC23_setFreq(hCodec, DSK6713_AIC23_FREQ_8KHZ);//set sample rate
	// Configure buffered serial ports for 32 bit operation
	// This allows transfer of both right and left channels in one read/write
	MCBSP_FSETS(SPCR1, RINTM, FRM);
	MCBSP_FSETS(SPCR1, XINTM, FRM);
	MCBSP_FSETS(RCR1, RWDLEN1, 32BIT);
	MCBSP_FSETS(XCR1, XWDLEN1, 32BIT);
	MCBSP_start(DSK6713_AIC23_DATAHANDLE, MCBSP_XMIT_START | MCBSP_RCV_START |
			MCBSP_SRGR_START | MCBSP_SRGR_FRAMESYNC, 220);// re-start data channel

	// initialize interrupt
	Uint32 CODECEventId;
	IRQ_globalDisable();		// disable interrupts
	CODECEventId=MCBSP_getXmtEventId(DSK6713_AIC23_codecdatahandle);// McBSP1 Xmit
	IRQ_setVecs(vectors);    	// point to the IRQ vector table
	IRQ_map(CODECEventId,11);	// Maps an event to a physical interrupt
	IRQ_reset(CODECEventId); 	// reset codec INT 11
	IRQ_globalEnable();			// Globally enables interrupts
	IRQ_nmiEnable();			// Enables the non-maskable interrupts
	IRQ_enable(CODECEventId);	// enable CODEC eventXmit INT11
	MCBSP_write(DSK6713_AIC23_DATAHANDLE,0);// start interrupt by outputting a sample
}
