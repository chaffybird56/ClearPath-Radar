/* Minimal stubs so RadarProcessor.c can compile on Linux CI (syntax check only). */
#ifndef CI_HOST_H
#define CI_HOST_H

#include <math.h>
#include <stdio.h>

typedef unsigned int Uint32;
typedef unsigned short Uint16;

#define interrupt
#define CHIP_6713

#define PTS 256
#define NUM_SAMPLES_PER_CHIRP 128
#define NUM_SAMPLES 256

typedef struct {
    float real;
    float imag;
} COMPLEX;

extern double dechirp_real[PTS];
extern double dechirp_img[PTS];
extern COMPLEX data_fft[PTS];
extern COMPLEX w_fft[PTS];
extern double data_fft_lib[PTS * 2];
extern double w_fft_lib[PTS * 2];

typedef void *TIMER_Handle;
typedef void *DSK6713_AIC23_CodecHandle;
typedef int DSK6713_AIC23_Config;

#define TIMER_DEVANY ((void *)0)
#define DSK6713_AIC23_DEFAULTCONFIG 0
#define DSK6713_AIC23_DATAHANDLE ((void *)0)
#define DSK6713_AIC23_FREQ_8KHZ 0
#define MCBSP_XMIT_START 0
#define MCBSP_RCV_START 0
#define MCBSP_SRGR_START 0
#define MCBSP_SRGR_FRAMESYNC 0

#define TIMER_configArgs(h, a, b, c) ((void)0)
#define TIMER_getCount(h) (0u)
#define TIMER_open(dev, f) ((TIMER_Handle)0)
#define TIMER_close(h) ((void)0)

#define DSK6713_init() ((void)0)
#define DSK6713_LED_init() ((void)0)
#define DSK6713_AIC23_openCodec(a, b) ((DSK6713_AIC23_CodecHandle)0)
#define DSK6713_AIC23_setFreq(h, f) ((void)0)
#define MCBSP_FSETS(a, b, c) ((void)0)
#define MCBSP_start(h, m, n) ((void)0)
#define MCBSP_getXmtEventId(h) (0u)
#define MCBSP_write(h, v) ((void)0)

#define IRQ_globalDisable() ((void)0)
#define IRQ_setVecs(v) ((void)0)
#define IRQ_map(e, n) ((void)0)
#define IRQ_reset(e) ((void)0)
#define IRQ_globalEnable() ((void)0)
#define IRQ_nmiEnable() ((void)0)
#define IRQ_enable(e) ((void)0)

void prepare_w_fft(COMPLEX *w, int N);
void prepare_w_fft_lib(double *w, int N);
void FFT(COMPLEX *data, COMPLEX *w, int N);
void DSPF_dp_cfftr2(int n, double *data, double *w, int flag);
void bit_rev_2(double *data, int n);

void vectors(void);

#endif
