#include "ci_host.h"

double dechirp_real[PTS];
double dechirp_img[PTS];
COMPLEX data_fft[PTS];
COMPLEX w_fft[PTS];
double data_fft_lib[PTS * 2];
double w_fft_lib[PTS * 2];

void prepare_w_fft(COMPLEX *w, int N) { (void)w; (void)N; }
void prepare_w_fft_lib(double *w, int N) { (void)w; (void)N; }
void FFT(COMPLEX *data, COMPLEX *w, int N) { (void)data; (void)w; (void)N; }
void DSPF_dp_cfftr2(int n, double *data, double *w, int flag)
{
    (void)n;
    (void)data;
    (void)w;
    (void)flag;
}
void bit_rev_2(double *data, int n) { (void)data; (void)n; }
void vectors(void) {}
