# Example measurement

`Sample_C1200_GITT_C-20_2.414mg.mpr`

One real GITT measurement from my master's thesis, included so the code can be tried out
without a cycler.

| | |
|---|---|
| Cell | hard carbon vs. sodium metal, half-cell (CR2032) |
| Sample | C1200 — coal-waste-derived hard carbon, carbonised at 1200 °C |
| Active mass | 2.414 mg (stored in the EC-Lab settings block, the app reads it automatically) |
| Electrode | 10 mm disc |
| Instrument | BioLogic VMP3, EC-Lab |
| Technique | GCPL used as GITT: 10-minute pulses at C/20 (1C = 250 mAh g⁻¹), 1-hour rests |
| Recorded | 2 October 2026 |
| Contents | 53 sodiation pulses, 95 h, 57 248 data points, 3.0 MB |
| Potential range | 0.005 – 2.484 V vs. Na⁺/Na |

This is the first sodiation half-cycle of a longer run; the file was cut after it so that it
stays small enough to keep in the repository.

## What the analysis gives

Pulse current 30.18 µA (C/20) · sodiation capacity 109.8 mAh g⁻¹ ·
*D*(Na⁺) median 5.3 × 10⁻¹⁰ cm² s⁻¹ with *L* = 20.5 µm (dense-carbon thickness).

The figures and the Excel table are in [`../example_output`](../example_output).

## It is not a perfect measurement

That is deliberate — the data checks are the point. The app reports:

1. 39 of the 52 used rests were still drifting faster than 3.6 mV h⁻¹ when the next pulse
   started, so those relaxed potentials are not equilibrium potentials.
2. The IR drop jumps from 41 mV to 71 mV between pulses 21 and 23 while the relaxed potential
   stays smooth, which points at the sodium counter electrode or a contact rather than at the
   hard carbon.

## Use

Measured by Bilal Dawood at TU Darmstadt. Shared under the same licence as the code, as a
demonstration dataset. If you use it for anything beyond trying the code out, please get in
touch first — the work it belongs to is not published yet.
