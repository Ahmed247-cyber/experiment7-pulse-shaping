# Speaker Script — Experiment 7: Pulse Shaping and the Nyquist Criterion

## Slide 1 — Experiment 7

Introduce the experiment as a comparison of pulse shapes and a verification of the Nyquist criterion. State that the key design question is how roll-off changes bandwidth and time-domain behavior.

## Slide 2 — Basic theory

Define the symbol period, normalized sinc function, raised-cosine pulse, and root-raised-cosine pulse. Explain that the RC response is the Nyquist response, while the RRC response is normally split between transmitter and receiver. Mention the removable singularities at the special time values and the explicit limiting values used in the implementation.

## Slide 3 — Nyquist criterion

Explain that zero intersymbol interference requires the overall pulse to equal one at the desired symbol time and zero at every other integer symbol interval. The RRC transmitter and matched RRC receiver cascade to an RC response, so the cascade is the response tested at symbol timing instants.

## Slide 4 — Experimental setup

Describe the 32-samples-per-symbol grid, the upsampled bipolar symbol impulse train, the deterministic 24-symbol stream, and the ±32-symbol filter span used for validation. Explain that the longer span is important for the alpha-zero sinc limit because sinc is theoretically infinite in duration.

## Slide 5 — Pulse comparison

Use the time-domain panel to compare compactness and sidelobe behavior. Use the frequency-domain panel to contrast the broad sidelobes of the rectangular pulse with the controlled transition of RC and RRC shaping. Note that finite FFTs introduce small spectral ripple.

## Slide 6 — Roll-off factor

Explain the expected effect before showing the results: alpha zero gives the narrowest bandwidth and slowest time decay; increasing alpha widens the transition band and accelerates time-domain decay. The tested values are 0, 0.25, 0.5, and 1.

## Slide 7 — Upsampling and pulse shaping

Point to the impulses spaced one symbol apart in the upper panel. Explain that convolution replaces isolated impulses with overlapping pulses. The overlap is intentional; the matched-filter cascade preserves the correct sample values at the decision instants.

## Slide 8 — RRC cascade

Explain that convolution in time corresponds to multiplication in frequency. Two matched RRC filters therefore produce the RC response. The plotted cascade shows unit gain at the origin and zeros at the integer symbol times.

## Slide 9 — Nyquist validation

Walk through the table. For each alpha, the k=0 sample is approximately one and the nonzero k samples are approximately zero. The alpha-zero case has the largest residual because truncating sinc creates the largest numerical error. The error decreases when the filter span is increased.

## Slide 10 — Bandwidth versus roll-off

State the theoretical relationship: B equals (1 plus alpha) divided by 2T. The bandwidth values are 0.5/T, 0.625/T, 0.75/T, and 1/T for the four tested roll-offs. Conclude that roll-off is a practical trade-off between spectral efficiency and time-domain localization.

## Slide 11 — Conclusions

Summarize three results: RC and cascaded RRC responses satisfy the Nyquist sampling condition; larger alpha increases bandwidth and reduces time-domain tails; finite truncation is the main source of residual error in the alpha-zero case.

## Slide 12 — Reproducibility

Identify the Python script, CSV validation tables, plots, and report. State that rerunning the script regenerates the numerical results and figures.
