# Roadside‑Clutter Radar — MTI Filtering on DSP 

> Suppress stationary roadside clutter near 0‑Hz Doppler so moving vehicles are emphasized. Compact FIR/IIR front‑end filters (MATLAB), two‑/three‑pulse MTI cancellers (DSP), and a windowed Doppler FFT complete the chain.

<!-- Scene/overview placeholder -->
<p align="center">
  <img width="640" height="280" alt="Scene overview" src="https://github.com/user-attachments/assets/4bb66543-2cd6-450a-b8f3-d4c7e62ca6f3" />
  <br/>
  <sub>Fig A — Roadside barrier scenario: clutter accumulates at 0 Hz; movers appear off zero.</sub>
</p>

---

## Rationale

Across sweeps, each range bin records a complex slow‑time sequence. Stationary reflectors retain phase and pile up at 0 Hz; moving reflectors accrue phase and shift energy away from 0 Hz. MTI performs sweep‑to‑sweep differencing to cancel near‑constant components before computing the Doppler FFT.

**Glossary.** MTI — Moving Target Indicator. FFT — Fast Fourier Transform. Range bin — one distance slice (one complex sample per sweep). Window — taper to control FFT sidelobes.

---

## Mathematical formulation

Fix a range bin and let $x[n]$ be the complex slow-time sample at sweep $n$.

### Two-pulse MTI (first difference)

$$
y[n] = x[n] - x[n-1]
$$

Frequency response $H_2(e^{j\omega}) = 1 - e^{-j\omega}$. Magnitude

$$
\lvert H_2(e^{j\omega})\rvert = 2\,\left\lvert\sin\left(\frac{\omega}{2}\right)\right\rvert,
$$

a high-pass with a notch at $\omega=0$.

### Three-pulse MTI (second difference)

$$
y[n] = x[n] - 2x[n-1] + x[n-2]
$$

Frequency response $H_3(e^{j\omega}) = (1 - e^{-j\omega})^2$. Magnitude

$$
\lvert H_3(e^{j\omega})\rvert = 4\,\sin^2\left(\frac{\omega}{2}\right),
$$

a wider, steeper notch around 0 Hz (more suppression of slowly varying clutter) with modest noise gain away from zero.
<!-- MTI response plot placeholder -->
<p align="center">
  <img width="647" height="363" alt="SCR-20250929-jtfb" src="https://github.com/user-attachments/assets/aed43adc-ae64-4b81-9ab3-4c939d6ce833" />
  <br/>
  <sub>Fig B — Example |H(e^{jω})| for two‑ and three‑pulse MTI.</sub>
</p>


### Windowed Doppler FFT

Hanning window

$$
w_h[n] = \tfrac{1}{2}\Big(1 - \cos\tfrac{2\pi n}{N-1}\Big),\qquad n=0,\dots,N-1
$$

Doppler spectrum

$$
X[k] = \sum_{n=0}^{N-1} \big(y[n]\,w_h[n]\big)\,e^{-j2\pi kn/N},\qquad k=0,\dots,N-1.
$$

The window reduces sidelobes so strong residuals near 0 Hz do not mask weak movers.

---

## Front‑end filter design (MATLAB)

A light FIR/IIR stage before MTI/FFT controls passband/stopband and dynamic range.

```matlab
% FIR (Dolph–Chebyshev) example
N   = 63;   As = 60;         % taps, stopband (dB)
fc  = 0.45;                   % normalized cutoff (0..1)
h   = fir1(N, fc, chebwin(N+1, As));

% IIR (Chebyshev II) example
Rs  = 60;  Wp = 0.45; Ws = 0.55;  % pass/stop edges
[n, Wn] = cheb2ord(Wp, Ws, 1, Rs);
[b, a]  = cheby2(n, Rs, Wn);

% Compare magnitude responses
[Hf, w] = freqz(h, 1, 4096);
[Hi, ~] = freqz(b, a, 4096);
plot(w/pi, 20*log10(abs(Hf)), 'LineWidth', 1.2); hold on;
plot(w/pi, 20*log10(abs(Hi)), 'LineWidth', 1.2); grid on;
xlabel('Normalized Frequency'); ylabel('Magnitude (dB)');
legend('FIR (Dolph–Cheb.)','IIR (Cheb. II)');
title('FIR vs IIR — Magnitude');
```
<!-- FIR/IIR plot placeholder -->
<p align="center">
  <img width="1267" height="407" alt="SCR-20250929-kdsl" src="https://github.com/user-attachments/assets/5474dd4c-50e6-41be-ba9c-c943922182b2" />
  <br/>
  <sub>Fig C — FIR vs IIR magnitude (illustrative).</sub>
</p>


```matlab
% Truncation vs windowing on FIR taps
N = 63; fc = 0.45;
h_box = fir1(N, fc, rectwin(N+1));
h_cheb = fir1(N, fc, chebwin(N+1, 60));
[Hb, w] = freqz(h_box, 1, 4096);
[Hc, ~] = freqz(h_cheb, 1, 4096);
plot(w/pi, 20*log10(abs(Hb)), 'LineWidth',1.2); hold on;
plot(w/pi, 20*log10(abs(Hc)), 'LineWidth',1.2); grid on;
legend('Rectangular truncation','Dolph–Cheb window');
title('FIR Truncation vs Windowing — Magnitude');
```
<!-- Trunc vs window plot placeholder -->
<p align="center">
  <img width="1229" height="400" alt="SCR-20250929-kdos" src="https://github.com/user-attachments/assets/37a0aa39-2288-4959-9d55-04d82b89a6a9" />
  <br/>
  <sub>Fig D — Truncation vs windowing (illustrative).</sub>
</p>


---

## DSP implementation (fixed‑point friendly)

Below are minimal C sources. Wire `fft_q15_inplace` to your platform library (DSPLib/CMSIS‑DSP/etc.) and select your Q‑format.

```c
// dsp/mti.h
#ifndef MTI_H
#define MTI_H
#include <stdint.h>
typedef int16_t q15_t; typedef int32_t q31_t;
typedef struct { q15_t x1, x2; } MTIState;
q15_t mti2(q15_t x, MTIState* s);   // y[n] = x[n] - x[n-1]
q15_t mti3(q15_t x, MTIState* s);   // y[n] = x[n] - 2x[n-1] + x[n-2]
#endif
```

```c
// dsp/mti.c
#include "mti.h"
static inline q15_t sat_q15(q31_t v){
    if(v >  32767) return  32767;
    if(v < -32768) return -32768;
    return (q15_t)v;
}
q15_t mti2(q15_t x, MTIState* s){
    q31_t y = (q31_t)x - (q31_t)s->x1;
    s->x1 = x;
    return sat_q15(y);
}
q15_t mti3(q15_t x, MTIState* s){
    q31_t y = (q31_t)x - ((q31_t)2*s->x1) + (q31_t)s->x2;
    s->x2 = s->x1; s->x1 = x;
    return sat_q15(y);
}
```

```c
// dsp/hann_table_q15.c
#include <stdint.h>
typedef int16_t q15_t; typedef int32_t q31_t;
/* Example: 256‑point Hanning in Q15. Replace with your N as needed. */
const q15_t HANN_Q15[256] = {
  0, 3, 13, 30, 53, 83, 119, 162, 211, 266, 328, 396, 470, 551, 637, 730,
  827, 930, 1038, 1151, 1268, 1390, 1516, 1646, 1780, 1917, 2058, 2202, 2349,
  /* ... fill to length as required ... */ 32764
};
```

```c
// dsp/fft_backend.c
#include <stdint.h>
/* Provide your platform FFT here; placeholder keeps linkage clean. */
void fft_q15_inplace(int16_t* buf, int N){
    (void)buf; (void)N; /* wire vendor DSPLib / CMSIS‑DSP / custom here */
}
```

**Processing order.**
1. Acquire \(x[n]\) per range bin.
2. Optional FIR/IIR front‑end.
3. MTI via `mti2` or `mti3` with per‑bin state.
4. Multiply by Hanning table, run \(N\)-point FFT, magnitude/threshold.

---

## Results

<!-- Scope BEFORE placeholder -->
<p align="center">
  <img width="639" height="361" alt="SCR-20250929-jtao" src="https://github.com/user-attachments/assets/4a977801-90fc-49de-9189-89ee6c348c0b" />
  <br/>
  <sub>Fig E — Doppler before MTI: strong 0‑Hz line from barrier.</sub>
</p>


<!-- Scope AFTER placeholder -->
<p align="center">
  <img width="647" height="363" alt="SCR-20250929-jtfb" src="https://github.com/user-attachments/assets/c8c4d261-7203-4b62-bd2c-7f80f0729354" />
  <br/>
  <sub>Fig F — Doppler after MTI + Hanning: clutter suppressed; movers visible at non‑zero bins.</sub>
</p>


---

## Tuning notes

- **Two‑ vs three‑pulse.** Two‑pulse maintains sensitivity to slow movers; three‑pulse deepens/widens the 0 Hz notch for harsher clutter.
- **Window choice.** Hanning is a solid default; heavier tapers lower sidelobes further at the cost of main‑lobe width.
- **FFT length \(N\).** Larger \(N\) improves velocity resolution but increases compute and coherent time.

---

## Repository layout

```
Roadside-Clutter-Radar/
├─ README.md
├─ dsp/
│  ├─ mti.h
│  ├─ mti.c
│  ├─ hann_table_q15.c
│  └─ fft_backend.c
└─ matlab/
   ├─ design_fir_vs_iir.m
   └─ trunc_vs_window.m
```

## License

MIT — see `LICENSE`.
