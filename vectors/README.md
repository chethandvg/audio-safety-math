# Golden vectors

Each `*.json` file here is a list of input/output cases that the Rust crate, the Python bindings and the WebAssembly bindings must all reproduce. The format is `golden-vectors@1` from [evidence-kernel](https://github.com/chethandvg/evidence-kernel).

The vectors are written by [`generate.py`](generate.py) from references that share no code with the crate:

- the decibel, level-sum, Leq and dose formulas evaluated directly in 50-digit arithmetic with mpmath (the crate factors sums around their largest term; the reference does not);
- `rms_dbfs` and `normalize` from the project's earlier toolkit, `lifesound.py`, copied with its own `box_lowpass` and `highpass` (not signal-core's, which the crate uses) and run with numpy, with two deliberate changes: `rms_dbfs` no longer adds 10⁻¹² before the logarithm, and `normalize` passes its sample rate to the high-pass; `peak_dbfs`, and `rms_dbfs` of samples whose squares underflow or overflow an `f64` (where numpy's version fails), in 50-digit arithmetic;
- the K-weighting coefficients built in 50-digit arithmetic from the 48 kHz coefficients published in ITU-R BS.1770 (tables 1 and 2), which are mapped back through the bilinear transform to two analogue prototypes and discretised again at each rate;
- the filter run with `scipy.signal.lfilter`;
- loudness measured with pyloudnorm 0.2.0's meter, its own filters replaced by those reference coefficients, and cross-checked against `bs1770`, a direct transcription of the BS.1770 formulas.

Nothing in the oracle calls the Rust code. Regenerate with:

```bash
python vectors/generate.py   # needs numpy, scipy, mpmath and pyloudnorm
```

| File | Functions | Tolerance |
|---|---|---|
| [`decibels.json`](decibels.json) | `amplitude_to_db`, `power_to_db`, `sum_incoherent_db`, `energetic_mean_db`, `leq` (1e-12); `db_to_amplitude`, `db_to_power` (1e-12 relative) | 1e-12 |
| [`dose.json`](dose.json) | `allowed_time`, `noise_dose`, `twa_from_dose` | 1e-12 |
| [`levels.json`](levels.json) | `rms_dbfs`, `peak_dbfs`, `normalize` | 1e-12 |
| [`filters.json`](filters.json) | `k_weighting` (1e-12; the published 48 kHz values 1e-13; 3364 Hz 1e-9 relative), `k_weight` (1e-12) | 1e-12 |
| [`loudness.json`](loudness.json) | `integrated_loudness`, `momentary_loudness`, `short_term_loudness` | 1e-9 LU |
| [`errors.json`](errors.json) | invalid input, in every function (the crate's own input rules, not oracle output) | — |

The six files hold 147 cases: 109 with expected outputs from the references above and 38 invalid inputs that must be rejected.

## Why these tolerances

Tolerances are `max(abs, rel · |expected|)`.

- **1e-12** (relative and absolute) for the decibel, dose and level functions: each is a few `log10`, `powf`, `exp2` or `sqrt` calls, each within about one unit in the last place (about 2 × 10⁻¹⁶ relative), on results of at most a few thousand dB; the remaining differences come from numpy's pairwise sums against the crate's running sums. The ratios returned by `db_to_amplitude` and `db_to_power` are checked relatively only, since an absolute 1e-12 would pass any value of 10⁻³⁰.
- **1e-13** for the published BS.1770 coefficients: they are given to 14 decimal places, and the crate's prototypes were derived from them, so the round trip is limited by that rounding (under 5 × 10⁻¹⁵) and the double-precision arithmetic of the bilinear transform.
- **1e-9 relative** for the K-weighting at 3364 Hz, the lowest rate: the shelf's corner sits just below the Nyquist frequency there, `tan(π·f₀/sr)` is about 41 900, and its condition number (about 6.6 × 10⁴) magnifies the rounding of `f₀/sr`.
- **1e-9 LU** for loudness: a block's mean square sums up to a few thousand squares (relative error at most about n·2⁻⁵³, under 10⁻¹² for these signals, so under 10⁻¹¹ LU), and the K-weighting's high-pass poles have radius 0.995 at 48 kHz (0.94 and 0.97 at the vectors' 4 and 8 kHz), which amplifies rounding in the recursion at most about two hundredfold (still under 10⁻¹³ relative per sample). 1e-9 LU leaves a margin of two orders of magnitude over both.

## What the cases cover

Every case names the claim in [`claims.yaml`](../claims.yaml) that it tests. Signals are flat arrays: multichannel signals are interleaved `[l0, r0, l1, r1, …]`, with the channel count as an argument (`normalize`, `k_weight`) or set by the list of channel weights (the loudness functions). To keep the files small, the loudness cases use low sample rates (4 and 8 kHz, above the K-weighting's minimum of 3364 Hz) and signals of 0.4 to 3.3 seconds rounded to four decimal places; they include signals whose blocks fall below the relative and the absolute gate, stereo, 5.0 and 5.1 with an LFE weight of 0, and a length that ends half-way through a 100 ms hop. Two references are not pyloudnorm, and say so in their notes: the 5.1 case (pyloudnorm's channel weights stop at five channels) and the short-term meter (pyloudnorm's 3 s blocks with an overlap of 1 − 1/30 start one sample early, because it truncates `3.0 · (j · (1 − overlap)) · sr` to a whole number), where `bs1770` agrees with pyloudnorm's integrated and momentary values within 10⁻¹¹ LU.

A few cases sit at the edges of the `f64` range: `leq` durations whose sum exceeds it, a `leq` segment at 6000 dB lasting 10⁻³³⁰ of the longest (a subnormal duration), a dose segment of 0 hours at 5000 dBA (which counts for nothing), 10⁻¹⁰ hours at 3160 dBA (whose `2^((L − L₀)/q)` overflows, though the dose does not), `rms_dbfs` of samples of 10⁻²⁰⁰ and 10³⁰⁰, and, among the invalid inputs, a constant that is not exact in binary (only rounding residue is left after `normalize`'s high-pass), a `k_weight` output and a loudness energy that overflow.

The 48 kHz conformance check against the EBU Tech 3341 test signals needs tens of seconds of audio per case, so it is not a vector file: `generate.py` writes the signal descriptions, the published expected values and pyloudnorm's readings to [`tests/ebu3341.json`](../tests/ebu3341.json), and [`tests/ebu3341.rs`](../tests/ebu3341.rs) checks the crate against them. Silence in the meters (−∞) and NaN cannot be written in JSON; [`tests/properties.rs`](../tests/properties.rs) and the binding tests check them.
