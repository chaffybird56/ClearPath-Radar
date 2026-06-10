# Roadside‑Clutter Radar — MTI Filtering on DSP 

[![CI](https://img.shields.io/github/actions/workflow/status/chaffybird56/ClearPath-Radar/ci.yml?label=CI&logo=github)](https://github.com/chaffybird56/ClearPath-Radar/actions/workflows/ci.yml)

**The problem in plain terms:** A radar looking down a road sees *everything* — guardrails, signs, parked cars, and the barrier along the shoulder. Those objects don't move, so their echoes look the same sweep after sweep. A car driving past *does* move, so its echo changes phase from sweep to sweep. This project shows how to **subtract the stationary stuff** (clutter) so **moving targets stand out** in the Doppler spectrum.

> **In one sentence:** Design and run a radar signal chain that cancels slow/stationary returns near 0 Hz (MTI), then FFTs what's left so movers appear at non‑zero Doppler bins — with MATLAB filter design up front and fixed‑point‑friendly C for the DSP.

<p align="center">
  <img width="640" height="280" alt="Scene overview" src="https://github.com/user-attachments/assets/4bb66543-2cd6-450a-b8f3-d4c7e62ca6f3" />
  <br/>
  <sub>Fig A — Roadside barrier scenario: clutter accumulates at 0 Hz; movers appear off zero.</sub>
</p>

## At a glance

| | |
|---|---|
| **Scenario** | FMCW / pulsed radar looking at a roadside barrier plus occasional moving vehicles. |
| **What clutter looks like** | Stationary reflectors → stable phase → energy piles up at **0 Hz Doppler**. |
| **What movers look like** | Changing range/phase across sweeps → energy shifts to **non‑zero Doppler bins**. |
| **What we do** | Optional FIR/IIR shaping → **MTI canceller** (sweep‑to‑sweep difference) → **Hanning window** → **Doppler FFT**. |
| **Where it runs** | MATLAB for design/plots; **TI TMS320C6713 DSK** for real‑time FFT + audio DAC playback; reference Q15 MTI C snippets for porting. |

## Signal chain (common language → technical)

```mermaid
flowchart LR
  A[Radar echo per range bin] --> B[Optional FIR/IIR]
  B --> C[MTI: subtract previous sweep]
  C --> D[Hanning window]
  D --> E[Doppler FFT]
  E --> F[Magnitude spectrum]
  F --> G[0 Hz clutter suppressed; movers visible]
```

1. **Acquire** — For each distance slice (*range bin*), record one complex sample per radar sweep (*slow time*).  
2. **Shape (optional)** — A light FIR or IIR filter sets passband/stopband before MTI. Designed in MATLAB.  
3. **MTI** — Compare this sweep to the last one (or last two). Anything that barely changed — the barrier — mostly cancels. Anything that moved — a car — survives.  
4. **Window + FFT** — Taper samples (Hanning) to control sidelobes, then FFT to get **velocity / Doppler** per bin.  
5. **Read the plot** — Strong line at 0 Hz before MTI = clutter. After MTI, that line drops and off‑zero peaks are movers.

**Glossary.** **MTI** — Moving Target Indicator. **FFT** — Fast Fourier Transform. **Range bin** — one distance slice (one complex sample per sweep). **Doppler** — how fast the target's phase changes across sweeps (related to radial velocity). **Window** — taper applied before FFT to reduce spectral leakage.

---

## Hardware setup

The on‑target demo targets the **Texas Instruments TMS320C6713 DSK** (Spectrum Digital DSK‑6713 EVM) — a floating‑point DSP starter kit common in university radar/DSP labs.

| Piece | Role in this project |
|-------|----------------------|
| **TMS320C6713** | 225 MHz C67x core; runs dechirp/FFT and streams magnitude bins in an interrupt service routine. |
| **AIC23 codec** (TLV320AIC23) | Audio DAC path: FFT magnitude per bin is scaled and played out at **8 kHz** so you can *hear* the spectrum evolve (left/right stereo on McBSP). |
| **Onboard USB emulator** | Debug/load via Code Composer Studio; target config in `targetConfigs/DSK6713.ccxml`. |
| **EMIF / SDRAM** | Linker map in `C6713.cmd` — program in internal IRAM; external SDRAM/flash/daughterboard regions mapped for DSK layout. |

**How `RadarProcessor.c` uses the board**

- **`initialize()`** — `DSK6713_init()` (PLL/EMIF), opens the AIC23 codec, sets **8 kHz** sample rate, configures **McBSP** for 32‑bit stereo frames, maps codec transmit to **CPU interrupt 11**.  
- **`main()`** — Loads precomputed dechirped chirp data (`data_dechirp`), prepares Hanning twiddles, runs a timed **FFT** (custom `FFT()` or `DSPF_dp_cfftr2` from DSPLib), prints cycle count, then idles.  
- **`c_int11()` ISR** — Each codec interrupt, outputs the magnitude of the next FFT bin through the DAC (with board‑specific sign inversion). Walking through bins at audio rate turns the range–Doppler slice into an audible sweep.

**Lab workflow (typical)**

1. Design / validate filters and MTI response in **MATLAB** (`Experiment.m`, README plots).  
2. Build and load `RadarProcessor.out` in **Code Composer Studio** with the DSK6713 target.  
3. Connect headphones or scope to the codec output; listen for the clutter notch and moving peaks after MTI.  
4. CI on GitHub only does a **host compile check** with `ci/stubs/` — not a full DSP link. See [docs/BUILD.md](docs/BUILD.md).

**Note:** Chirp/dechirp sample arrays ship as generated headers (`data_dechirp.h`, etc.) from the course radar lab pipeline — the repo demonstrates the **processing chain and MTI concept**, not a full RF front‑end schematic.

---

## Rationale (technical)

Across sweeps, each range bin records a complex slow‑time sequence. Stationary reflectors retain phase and pile up at 0 Hz; moving reflectors accrue phase and shift energy away from 0 Hz. MTI performs sweep‑to‑sweep differencing to cancel near‑constant components before computing the Doppler FFT.

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

The window reduces sidelobes so strong residuals near 0 Hz do not mask weak movers.

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

<p align="center">
  <img width="1229" height="400" alt="SCR-20250929-kdos" src="https://github.com/user-attachments/assets/37a0aa39-2288-4959-9d55-04d82b89a6a9" />
  <br/>
  <sub>Fig D — Truncation vs windowing (illustrative).</sub>
</p>


---

## DSP implementation (fixed‑point friendly)

Reference C for porting MTI to Q15 on any DSP (wire `fft_q15_inplace` to DSPLib/CMSIS‑DSP/etc.).

```c
// mti.h (reference)
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
// mti.c (reference)
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

**Processing order (per range bin).**
1. Acquire $x[n]$ each sweep.  
2. Optional FIR/IIR front‑end.  
3. MTI via `mti2` or `mti3` with per‑bin state.  
4. Multiply by Hanning table, run $N$‑point FFT, magnitude/threshold.

On the **C6713 DSK**, `RadarProcessor.c` implements steps 4+ playback using floating‑point FFT and preloaded dechirped data; MTI can be inserted in the slow‑time loop before the FFT stage when extending the lab.

---

## Results

<p align="center">
  <img width="639" height="361" alt="SCR-20250929-jtao" src="https://github.com/user-attachments/assets/4a977801-90fc-49de-9189-89ee6c348c0b" />
  <br/>
  <sub>Fig E — Doppler before MTI: strong 0 Hz line from barrier.</sub>
</p>

<p align="center">
  <img width="647" height="363" alt="SCR-20250929-jtfb" src="https://github.com/user-attachments/assets/c8c4d261-7203-4b62-bd2c-7f80f0729354" />
  <br/>
  <sub>Fig F — Doppler after MTI + Hanning: clutter suppressed; movers visible at non‑zero bins.</sub>
</p>


---

## Tuning notes

- **Two‑ vs three‑pulse.** Two‑pulse maintains sensitivity to slow movers; three‑pulse deepens/widens the 0 Hz notch for harsher clutter.  
- **Window choice.** Hanning is a solid default; heavier tapers lower sidelobes further at the cost of main‑lobe width.  
- **FFT length $N$.** Larger $N$ improves velocity resolution but increases compute and coherent integration time.  
- **On hardware.** If audio playback sounds saturated, rescale in `convert_double_to_short_for_output()`; codec output is inverted on this board (`sig_out = -1*sig_out`).

---

## Repository layout

```
ClearPath-Radar/
├─ README.md
├─ RadarProcessor.c      # C6713 DSK: FFT + AIC23 DAC playback (ISR)
├─ Experiment.m          # MATLAB DFT / filter design experiments
├─ C6713.cmd             # Linker command file (IRAM, SDRAM, EMIF)
├─ Vectors_intr.asm      # Interrupt vector table
├─ targetConfigs/        # CCS target: DSK6713 + onboard USB emulator
├─ docs/BUILD.md         # Build notes (CCS vs CI stubs)
├─ ci/stubs/             # Host-side headers for GitHub Actions compile check
└─ Debug/                # Code Composer Studio build artifacts
```

## License

MIT — see `LICENSE`.
