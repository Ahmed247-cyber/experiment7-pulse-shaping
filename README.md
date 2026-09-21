# Experiment 7 — Pulse Shaping and the Nyquist Criterion

This repository contains a reproducible raised-cosine/root-raised-cosine pulse-shaping experiment.

## Contents

- `experiment7_pulse_shaping.py` — Python simulation and plotting script.
- `experiment7_report.md` — theory, results, observations, and Nyquist validation.
- `nyquist_zero_crossings.csv` — integer-symbol samples of the cascaded RRC response.
- `bandwidth_rolloff.csv` — theoretical bandwidth and numerical diagnostics.
- `01_pulse_comparison.png` through `05_bandwidth_vs_rolloff.png` — required visualizations.
## Run

```bash
python3 experiment7_pulse_shaping.py
```

The script requires NumPy, SciPy, and Matplotlib.
