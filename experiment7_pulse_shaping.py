#!/usr/bin/env python3
"""Experiment 7: Pulse shaping and the Nyquist criterion.

Generates mathematical RC/RRC filters, time/frequency plots, a pulse-shaped
symbol stream, bandwidth-vs-rolloff data, and the mandatory integer-symbol
sample table for the cascaded RRC response.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import freqz

OUT = Path('/home/ubuntu/experiment7_output')
OUT.mkdir(exist_ok=True)
SPS = 32                 # samples per symbol
SPAN = 32                # half-span in symbols; larger span tests alpha=0 fairly
N = 2 * SPAN * SPS + 1
T = 1.0
ALPHAS = [0.0, 0.25, 0.5, 1.0]

plt.rcParams.update({'figure.dpi': 150, 'savefig.dpi': 180, 'font.size': 10,
                     'axes.grid': True, 'grid.alpha': 0.25})

def time_axis():
    return np.arange(-SPAN*SPS, SPAN*SPS + 1) / SPS

def rc_pulse(t, alpha, T=1.0):
    """Raised-cosine impulse response, with removable singularities explicit."""
    x = t / T
    if alpha == 0:
        return np.sinc(x)
    with np.errstate(divide='ignore', invalid='ignore'):
        y = np.sinc(x) * np.cos(np.pi * alpha * x) / (1 - (2 * alpha * x)**2)
    singular = np.isclose(np.abs(x), 1/(2*alpha), atol=1e-12)
    # lim h(t) at t=+-T/(2 alpha) = (alpha/2) sin(pi/(2 alpha))
    y[singular] = (alpha/2) * np.sin(np.pi/(2*alpha))
    return y

def rrc_pulse(t, alpha, T=1.0):
    """Root-raised-cosine impulse response with t=0 and t=+-T/(4a) limits."""
    x = t / T
    if alpha == 0:
        return np.sinc(x)
    with np.errstate(divide='ignore', invalid='ignore'):
        y = (np.sin(np.pi*x*(1-alpha)) + 4*alpha*x*np.cos(np.pi*x*(1+alpha))) / (np.pi*x*(1-(4*alpha*x)**2))
    at0 = np.isclose(x, 0, atol=1e-12)
    y[at0] = 1 - alpha + 4*alpha/np.pi
    singular = np.isclose(np.abs(x), 1/(4*alpha), atol=1e-12)
    # removable limit at x=+-1/(4 alpha)
    y[singular] = (alpha/np.sqrt(2)) * ((1+2/np.pi)*np.sin(np.pi/(4*alpha)) + (1-2/np.pi)*np.cos(np.pi/(4*alpha)))
    return y

def rectangular_pulse(t, T=1.0):
    return (np.abs(t) < T/2).astype(float)

def spectrum(h, nfft=32768):
    H = np.fft.fftshift(np.fft.fft(np.fft.ifftshift(h), nfft)) / SPS
    f = np.fft.fftshift(np.fft.fftfreq(nfft, d=1/SPS))
    mag = 20*np.log10(np.maximum(np.abs(H)/np.max(np.abs(H)), 1e-8))
    return f, mag

# Base pulses comparison.
t = time_axis()
pulses = {
    'Rectangular': rectangular_pulse(t),
    'Sinc': np.sinc(t),
    'Raised cosine (a=0.5)': rc_pulse(t, .5),
    'Root-raised cosine (a=0.5)': rrc_pulse(t, .5),
}
fig, ax = plt.subplots(2, 1, figsize=(9, 7), constrained_layout=True)
for name, h in pulses.items():
    ax[0].plot(t, h, label=name)
    f, m = spectrum(h)
    keep = np.abs(f) <= 2.0
    ax[1].plot(f[keep], m[keep], label=name)
ax[0].set(xlim=(-4, 4), ylabel='Amplitude', title='Pulse impulse responses')
ax[1].set(xlim=(-2, 2), ylim=(-80, 3), xlabel='Frequency (cycles/symbol)', ylabel='Magnitude (dB)', title='Pulse frequency responses')
ax[0].legend(ncol=2); ax[1].legend(ncol=2)
fig.savefig(OUT/'01_pulse_comparison.png'); plt.close(fig)

# RC roll-off comparison: expected effect is documented in report.
fig, ax = plt.subplots(2, 1, figsize=(9, 7), constrained_layout=True)
for a in ALPHAS:
    h = rc_pulse(t, a)
    ax[0].plot(t, h, label=f'alpha={a}')
    f, m = spectrum(h)
    keep = np.abs(f) <= 1.25
    ax[1].plot(f[keep], m[keep], label=f'alpha={a}')
ax[0].set(xlim=(-4, 4), ylabel='Amplitude', title='Raised-cosine impulse response vs roll-off')
ax[1].set(xlim=(-.8, .8), ylim=(-80, 3), xlabel='Frequency (cycles/symbol)', ylabel='Magnitude (dB)', title='Raised-cosine frequency response vs roll-off')
ax[0].legend(ncol=2); ax[1].legend(ncol=2)
fig.savefig(OUT/'02_rc_rolloff.png'); plt.close(fig)

# Deterministic upsampled symbol impulse train and pulse-shaped stream.
rng = np.random.default_rng(7)
symbols = rng.choice([-1.0, 1.0], size=24)
upsampled = np.zeros(len(symbols)*SPS)
upsampled[::SPS] = symbols
h_rrc = rrc_pulse(t, .5)
stream = np.convolve(upsampled, h_rrc, mode='same')
fig, ax = plt.subplots(2, 1, figsize=(10, 6), constrained_layout=True)
ax[0].stem(np.arange(len(upsampled))/SPS, upsampled, linefmt='C0-', markerfmt='C0.', basefmt='k-', label='upsampled impulse train')
ax[0].set(xlim=(0, 24), ylabel='Symbol amplitude', title='Upsampled symbol impulse train')
ax[1].plot(np.arange(len(stream))/SPS, stream, label='RRC pulse-shaped stream')
ax[1].set(xlim=(0, 24), xlabel='Time (symbols)', ylabel='Amplitude', title='Pulse-shaped stream (alpha=0.5)')
for a in ax: a.legend()
fig.savefig(OUT/'03_pulse_shaped_stream.png'); plt.close(fig)

# Cascade transmitter and receiver RRC filters. In continuous theory this is RC;
# finite truncation and discrete convolution are retained and quantified.
fig, ax = plt.subplots(2, 1, figsize=(9, 7), constrained_layout=True)
for a in ALPHAS:
    h = rrc_pulse(t, a)
    cascade = np.convolve(h, h) / SPS  # discrete approximation to continuous convolution
    tc = np.arange(cascade.size)/SPS - 2*SPAN
    ax[0].plot(tc, cascade, label=f'alpha={a}')
    f, m = spectrum(cascade)
    keep = np.abs(f) <= 1.25
    ax[1].plot(f[keep], m[keep], label=f'alpha={a}')
ax[0].set(xlim=(-4, 4), ylabel='Amplitude', title='Cascaded transmitter + receiver RRC response (RC)')
ax[1].set(xlim=(-.8, .8), ylim=(-80, 3), xlabel='Frequency (cycles/symbol)', ylabel='Magnitude (dB)', title='Cascaded RRC frequency response')
ax[0].legend(ncol=2); ax[1].legend(ncol=2)
fig.savefig(OUT/'04_rrc_cascade.png'); plt.close(fig)

# Mandatory validation table: sample overall RC response at integer symbol intervals.
rows = []
for a in ALPHAS:
    h = rrc_pulse(t, a)
    cascade = np.convolve(h, h) / SPS  # discrete approximation to continuous convolution
    tc = np.arange(cascade.size)/SPS - 2*SPAN
    samples = []
    for k in range(-6, 7):
        idx = np.argmin(np.abs(tc-k))
        samples.append(float(cascade[idx]))
    rows.append((a, samples))
with open(OUT/'nyquist_zero_crossings.csv', 'w', newline='') as fp:
    w = csv.writer(fp)
    w.writerow(['alpha'] + [f'k={k}' for k in range(-6,7)])
    for a, samples in rows: w.writerow([a] + [f'{v:.10e}' for v in samples])

# Numerical diagnostics and bandwidth-vs-rolloff. Theoretical one-sided Nyquist
# bandwidth is B=(1+alpha)/(2T); report measured -3 dB occupied edge too.
metrics = []
for a in ALPHAS:
    h = rc_pulse(t, a)
    f, mag = spectrum(h)
    pos = f >= 0
    fp, mp = f[pos], mag[pos]
    # Edge where response first drops below -60 dB after the passband.
    idx = np.where((fp > .45) & (mp < -60))[0]
    edge = float(fp[idx[0]]) if len(idx) else np.nan
    max_off = max(abs(v) for _, samples in rows if _ == a for k,v in zip(range(-6,7), samples) if k != 0)
    metrics.append((a, (1+a)/(2*T), edge, max_off))
with open(OUT/'bandwidth_rolloff.csv', 'w', newline='') as fp:
    w = csv.writer(fp); w.writerow(['alpha','theoretical_nyquist_bandwidth_cycles_per_symbol','measured_approx_edge_at_-60dB','max_abs_integer_symbol_off_center_sample'])
    for row in metrics: w.writerow([f'{v:.10e}' if isinstance(v,float) else v for v in row])
fig, ax = plt.subplots(figsize=(8, 4.8), constrained_layout=True)
a = np.array([r[0] for r in metrics]); bw = np.array([r[1] for r in metrics])
ax.plot(a, bw, 'o-', label='Theory: (1+alpha)/(2T)')
ax.set(xlabel='Roll-off factor alpha', ylabel='One-sided bandwidth (cycles/symbol)', title='Bandwidth versus roll-off factor')
ax.set_xticks(ALPHAS); ax.legend()
fig.savefig(OUT/'05_bandwidth_vs_rolloff.png'); plt.close(fig)

# Text summary used by the report.
with open(OUT/'metrics.txt', 'w') as fp:
    fp.write('Experiment 7 numerical diagnostics\n')
    fp.write(f'SPS={SPS}, truncated pulse span=+/-{SPAN} symbols\n\n')
    for a, bw, edge, off in metrics:
        fp.write(f'alpha={a}: theory bandwidth={bw:.6f}, approx -60 dB edge={edge:.6f}, max |off-center integer sample|={off:.3e}\n')
print(f'Wrote figures and CSV tables to {OUT}')
