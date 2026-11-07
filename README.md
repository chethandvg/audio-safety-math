# Audio safety math

**Decibel arithmetic, noise-dose rules and loudness metering in Rust, Python and WebAssembly, with Lean 4 proofs of the dose rule.**

[![CI](https://github.com/chethandvg/audio-safety-math/actions/workflows/ci.yml/badge.svg)](https://github.com/chethandvg/audio-safety-math/actions/workflows/ci.yml)
[![Proofs](https://github.com/chethandvg/audio-safety-math/actions/workflows/proofs.yml/badge.svg)](https://github.com/chethandvg/audio-safety-math/actions/workflows/proofs.yml)
[![Licence](https://img.shields.io/badge/licence-Apache--2.0-blue)](LICENSE)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24044023.svg)](https://doi.org/10.5281/zenodo.24044023)

Part of **[Frequencies of Life](https://chethandvg.github.io/frequencies-of-life/)**, open code for the rhythms of life in which every claim says how it is known: proved, verified, tested, referenced or assumed.

> **Not medical advice.** This library computes decibels, noise doses and loudness from formulas and published reference values. It does not measure sound at anyone's ear, it is not a hearing test, and a dose below 100 % does not mean a sound is safe for a particular person. For hearing concerns, ask an audiologist or a doctor; for workplace noise, follow the rules that apply where you work.

## Claims and evidence

<!-- claims-table:start -->
| Claim | Level | Statement | Evidence |
|---|---|---|---|
| `safety:DB-PRODUCT` | P | For amplitude ratios a, b > 0, 20·log₁₀(a·b) = 20·log₁₀ a + 20·log₁₀ b: gains in a chain multiply, so their decibels add. | [`AudioSafetyMath.amplitudeDb_mul`](proofs/AudioSafetyMath/Decibel.lean) |
| `safety:INCOHERENT-SUM-BOUNDS` | P | For n ≥ 1 levels Lᵢ, max Lᵢ ≤ 10·log₁₀ Σ 10^(Lᵢ/10) ≤ max Lᵢ + 10·log₁₀ n: unrelated sources together are at least as loud as the loudest one and at most 10·log₁₀ n dB louder. | [`AudioSafetyMath.sumDb_bounds`](proofs/AudioSafetyMath/Decibel.lean) |
| `safety:EXCHANGE-RATE` | P | With allowed time T(L) = T_c / 2^((L − L₀)/q) for criterion level L₀, exchange rate q ≠ 0 and criterion duration T_c, T(L + q) = T(L)/2, so a dose of t hours at L equals a dose of t/2 hours at L + q. | [`AudioSafetyMath.dose_half_time_add_exchange`](proofs/AudioSafetyMath/Dose.lean) |
| `safety:THREE-DB-IS-ABOUT-DOUBLING` | P | 3.0102 < 10·log₁₀ 2 < 3.0103: doubling the power adds just over 3 dB, so a 3 dB exchange rate approximates the equal-energy rule. | [`AudioSafetyMath.three_db_about_doubling`](proofs/AudioSafetyMath/Decibel.lean) |
| `safety:TWA-INVERTS-DOSE` | P | A whole criterion duration T_c ≠ 0 at a steady level L is a dose of 100·2^((L − L₀)/q) percent, and the time-weighted average L₀ + q·log₂(D/100) of that dose is L again (q ≠ 0). | [`AudioSafetyMath.twa_dose`](proofs/AudioSafetyMath/Dose.lean) |
| `safety:NIOSH-REL` | R | The default dose parameters are the NIOSH recommended exposure limit for occupational noise: 85 dBA as an 8-hour time-weighted average, with a 3 dB exchange rate (every 3 dB more halves the allowed time). | NIOSH 1998, Criteria for a recommended standard: occupational noise exposure, revised criteria 1998, DHHS (NIOSH) publication 98-126, Chapter 1, Section 1.1 ([doi:10.26616/NIOSHPUB98126](https://doi.org/10.26616/NIOSHPUB98126)) |
| `safety:DB-CONVERSIONS` | T | amplitude_to_db, db_to_amplitude, power_to_db and db_to_power match the formulas evaluated in 50-digit arithmetic within 1e-12, decibels of products add within 1e-12 dB in the code, conversions round-trip within 1e-13 relative, and 10·log₁₀ 2 in the code lies between 3.0102 and 3.0103. | [`vectors/decibels.json`](vectors/decibels.json), [`tests/properties.rs`](tests/properties.rs), [`src/db.rs`](src/db.rs) |
| `safety:LEVEL-SUMS` | T | sum_incoherent_db, energetic_mean_db and leq match the formulas evaluated in 50-digit arithmetic within 1e-12; in the floating-point code the incoherent sum always lies within the proved bounds, and levels from −5000 to 7000 dB give finite results, as does leq of a segment at 6000 dB lasting 10⁻³³⁰ of the longest (which matches 50-digit arithmetic), because each sum is factored around its largest term, with the durations entering as logarithms. | [`vectors/decibels.json`](vectors/decibels.json), [`tests/properties.rs`](tests/properties.rs), [`src/db.rs`](src/db.rs) |
| `safety:LEVELS-MATCH-REFERENCE` | T | rms_dbfs and normalize (on signal-core's zero-phase 12 Hz high-pass) match lifesound.py's rms_dbfs and normalize, run in numpy with lifesound.py's own box_lowpass and highpass, within 1e-12, and peak_dbfs, and rms_dbfs of samples whose squares underflow or overflow an f64, match 20·log₁₀ max\|x\| and 20·log₁₀ √(mean(x²)) in 50-digit arithmetic; rms_dbfs no longer adds 1e-12 before the logarithm, so silence is an error rather than −240 dBFS. normalize reaches the RMS target unless the peak ceiling binds, and then puts the peak at the ceiling; a signal left with only rounding residue after the high-pass (a constant, exact in binary or not) is an error rather than noise scaled to the target, and its gain overflows only when the result is too large for an f64 (not for a ceiling that does not bind, nor to lift subnormal samples); samples so large that the high-pass's running sums overflow, which takes a peak of at least about 1.8 × 10³⁰⁸/(len(x) + L) for the box length L, are an error on x, even when the result would fit. | [`vectors/levels.json`](vectors/levels.json), [`tests/properties.rs`](tests/properties.rs), [`src/db.rs`](src/db.rs) |
| `safety:DOSE-MATCHES-REFERENCE` | T | allowed_time, noise_dose and twa_from_dose match the formulas evaluated in 50-digit arithmetic within 1e-12 (with the NIOSH parameters, 88 dBA is allowed exactly 4 h and 100 dBA 15 min); in the floating-point code one exchange rate more halves the allowed time and keeps the dose of half the time within 1e-13 relative, and the TWA of a full day at a steady level returns that level within 1e-12 dB; the shares of the dose are summed as powers of two around the largest, so 10⁻¹⁰ hours at 3160 dBA, whose 2^((L − L₀)/q) overflows an f64, still gives its dose of 4.49 × 10²⁹⁹ %. | [`vectors/dose.json`](vectors/dose.json), [`tests/properties.rs`](tests/properties.rs), [`src/dose.rs`](src/dose.rs) |
| `safety:K-WEIGHTING-COEFFICIENTS` | T | k_weighting builds both K-weighting stages from fixed analogue prototypes by the bilinear transform prewarped at each stage's corner frequency; at 48 kHz it reproduces the coefficients published in ITU-R BS.1770 within 1e-14, at other rates it matches the same construction in 50-digit arithmetic within 1e-12 (1e-9 at the lowest rate, 3364 Hz), and below 6 kHz its response at 44.1, 96 and 192 kHz stays within 0.01 dB of the 48 kHz filter. | [`vectors/filters.json`](vectors/filters.json), [`tests/properties.rs`](tests/properties.rs), [`src/loudness.rs`](src/loudness.rs) |
| `safety:K-WEIGHT-FILTER` | T | k_weight filters each channel of an interleaved signal through both stages from rest and matches scipy.signal.lfilter with the reference coefficients within 1e-12. | [`vectors/filters.json`](vectors/filters.json) |
| `safety:LOUDNESS-MATCHES-PYLOUDNORM` | T | integrated_loudness and momentary_loudness match pyloudnorm 0.2.0's meter run with the reference K-weighting within 1e-9 LU, for mono, stereo and 5.0 signals and signals that exercise both gates; short_term_loudness and a 5.1 case with an LFE weight of 0 match a direct transcription of the BS.1770 formulas within 1e-9 LU; a gain of g adds 20·log₁₀ g LU and a copied channel adds 10·log₁₀ 2 within 1e-10. | [`vectors/loudness.json`](vectors/loudness.json), [`tests/properties.rs`](tests/properties.rs), [`src/loudness.rs`](src/loudness.rs) |
| `safety:EBU-3341` | T | On the EBU Tech 3341 test signals 1 to 6 and 9 (1 kHz sines at 48 kHz), integrated_loudness reads −22.993 (case 1), −32.993 (2), −23.014 (3 and 4), −22.979 (5) and −23.016 (6) LUFS, momentary and short-term loudness stay within 0.007 LU of −23 or −33 LUFS (cases 1 and 2) and short-term within 0.014 LU of −23 LUFS (case 9): all inside the published ±0.1 LU, and the integrated readings within 1e-6 LU of pyloudnorm's. | [`tests/ebu3341.rs`](tests/ebu3341.rs), [`tests/ebu3341.json`](tests/ebu3341.json) |
| `safety:INPUT-VALIDATION` | T | Every function rejects NaN and infinite input, silence where a single level is asked for (a zero ratio, an all-zero signal, a signal left with only rounding residue after normalize's high-pass, every loudness block below the absolute gate), negative durations and weights, mismatched lists, buffers that are not whole frames, sample rates below 3364 Hz for loudness, and results that overflow an f64, including a filtered signal or a loudness energy that overflows on the way; Rust returns SafetyError, Python raises ValueError (OverflowError for a negative or too-large integer argument such as a sample rate, TypeError for a non-integer one) and JavaScript throws an Error. The meters return −∞ for a silent window. | [`vectors/errors.json`](vectors/errors.json), [`tests/properties.rs`](tests/properties.rs), [`src/lib.rs`](src/lib.rs), [`src/db.rs`](src/db.rs), [`src/dose.rs`](src/dose.rs), [`src/loudness.rs`](src/loudness.rs) |
| `safety:IMPLEMENTATIONS-AGREE` | T | The Rust crate, the Python bindings and the WebAssembly bindings all reproduce every golden vector. Every expected output comes from an independent reference (mpmath, lifesound.py with numpy, scipy, pyloudnorm, or a transcription of BS.1770), never from the Rust output; the error cases in errors.json are the crate's own input rules. | [`vectors/generate.py`](vectors/generate.py), [`tests/vectors.rs`](tests/vectors.rs), [`bindings/python/tests/test_vectors.py`](bindings/python/tests/test_vectors.py), [`bindings/wasm/tests/vectors.test.mjs`](bindings/wasm/tests/vectors.test.mjs) |
<!-- claims-table:end -->

Levels: **P** machine-checked proof (Lean 4) · **V** verified numerics · **T** tested · **R** referenced (DOI or PMID) · **M** model assumption.

## Use

Three groups of functions, with the same names and positional arguments in Rust, Python and JavaScript:

- **Decibels:** `amplitude_to_db` and `db_to_amplitude` (20·log₁₀), `power_to_db` and `db_to_power` (10·log₁₀), `sum_incoherent_db` (unrelated sources together: their powers add), `energetic_mean_db` and `leq` (the equal-energy average of levels, the latter for timed segments), and the signal levels `rms_dbfs`, `peak_dbfs` and `normalize`, ported from the project's earlier toolkit `lifesound.py`.
- **Noise dose:** `allowed_time(level_dba, criterion_db, exchange_db, criterion_hours)` = `criterion_hours / 2^((level − criterion)/exchange)`, the daily `noise_dose` in percent (`100 · Σ Cᵢ/Tᵢ`), and `twa_from_dose`, the steady level that gives the same dose. The usual parameters are NIOSH's recommended exposure limit, 85 dBA as an 8-hour time-weighted average with a 3 dB exchange rate (see [Exposure limits referenced](#exposure-limits-referenced); exported as `NIOSH_CRITERION_DB`, `NIOSH_EXCHANGE_DB` and `CRITERION_HOURS` in Rust and Python, and the defaults in Python and JavaScript); OSHA's 90 dBA with a 5 dB exchange rate is the same formula with other parameters.
- **Loudness** per ITU-R BS.1770: `k_weighting(sr)` (the two filter stages), `k_weight` (a signal through them), the gated `integrated_loudness` in LUFS, and `momentary_loudness` (400 ms windows) and `short_term_loudness` (3 s windows), one value every 100 ms. Signals are interleaved frames (`[l0, r0, l1, r1, …]`) and the list of channel weights sets the channel count: `[1, 1]` for stereo, `[1, 1, 1, 1.41, 1.41]` for 5.0 (L, R, C, Ls, Rs), 0 for an LFE channel.

**Rust**

```toml
[dependencies]
audio-safety-math = { git = "https://github.com/chethandvg/audio-safety-math", tag = "v0.1.0" }
signal-core = { git = "https://github.com/chethandvg/signal-core", tag = "v0.1.0" } # only for osc in this example
```

```rust
use audio_safety_math::{
    allowed_time, integrated_loudness, leq, noise_dose, sum_incoherent_db, twa_from_dose,
};
use signal_core::osc;

// Two machines at 85 dBA each are 88 dBA together: 4 hours a day instead of 8.
let both = sum_incoherent_db(&[85.0, 85.0])?; // 88.0103 dB
let hours = allowed_time(both, 85.0, 3.0, 8.0)?; // 3.99 h (NIOSH: 85 dBA, 3 dB, 8 h)

// A working day of 2 h at 91 dBA and 6 h at 80 dBA.
let dose = noise_dose(&[91.0, 80.0], &[2.0, 6.0], 85.0, 3.0, 8.0)?; // 123.6 %
let twa = twa_from_dose(dose, 85.0, 3.0)?; // 85.92 dBA
let level = leq(&[91.0, 80.0], &[2.0, 6.0])?; // 85.91 dBA, the equal-energy average

// One second of a 1 kHz sine at -23 dBFS in both channels of a stereo file.
let tone = osc(1000.0, 1.0, 48_000, 0.0)?;
let stereo: Vec<f64> = tone.iter().flat_map(|&v| [0.0708 * v; 2]).collect();
let lufs = integrated_loudness(&stereo, &[1.0, 1.0], 48_000)?; // -22.99 LUFS
```

**Python** (builds the Rust core with maturin; signals come back as numpy arrays)

```bash
pip install "git+https://github.com/chethandvg/audio-safety-math@v0.1.0#subdirectory=bindings/python"
```

```python
import audio_safety_math as asm
import numpy as np

asm.sum_incoherent_db([85.0, 85.0])  # 88.0103: two equal sources are 3.01 dB louder
asm.allowed_time(88.0)  # 4.0 hours with the NIOSH defaults (85 dBA, 3 dB, 8 h)
asm.allowed_time(95.0, criterion_db=90.0, exchange_db=5.0)  # 4.0 hours by OSHA's 90 dBA / 5 dB rule
dose = asm.noise_dose([91.0, 80.0], [2.0, 6.0])  # 123.6 (percent): 2 h at 91 dBA, 6 h at 80 dBA
asm.twa_from_dose(dose)  # 85.92 dBA

sr = 48000
tone = 10 ** (-23 / 20) * np.sin(2 * np.pi * 1000 * np.arange(sr) / sr)  # 1 s, 1 kHz, -23 dBFS
stereo = np.repeat(tone, 2)  # interleaved [l0, r0, l1, r1, ...], the same in both channels
asm.integrated_loudness(stereo, [1.0, 1.0], sr)  # -22.99 LUFS
asm.momentary_loudness(stereo, [1.0, 1.0], sr)  # numpy array of 7 values, one per 100 ms
asm.k_weighting(48000)  # [[1.53512485958697, ...], [1.0, -2.0, 1.0, ...]]: BS.1770's coefficients
asm.rms_dbfs([0.0, 0.0])  # ValueError: x must be a signal that is not silent: ...
```

**JavaScript / WebAssembly** (build steps in [bindings/wasm](bindings/wasm))

```js
import init, { allowed_time, noise_dose, integrated_loudness } from "./pkg/web/audio_safety_math_wasm.js";

await init();
allowed_time(88);                  // 4 (hours; omitted arguments take the NIOSH rule)
allowed_time(95, 90, 5);           // 4 (OSHA's parameters)
noise_dose([91, 80], [2, 6]);      // 123.62 (percent)
const sr = 48000;
const tone = Float64Array.from({ length: sr }, (_, i) => 10 ** (-23 / 20) * Math.sin((2 * Math.PI * 1000 * i) / sr));
const stereo = Float64Array.from({ length: 2 * sr }, (_, i) => tone[i >> 1]);
integrated_loudness(stereo, [1, 1], sr);  // -22.99 (LUFS)
```

In Python the dose functions take keyword defaults (`criterion_db=85.0`, `exchange_db=3.0`, `criterion_hours=8.0`) and the functions that take a signal take none; in JavaScript a trailing dose argument left `undefined` takes the NIOSH value, and channel counts and sample rates are numbers that must be whole. `k_weighting` returns two lists of `[b0, b1, b2, a1, a2]` (`a0 = 1`), as two `Float64Array`s in JavaScript. The constants (`NIOSH_CRITERION_DB`, `NIOSH_EXCHANGE_DB`, `CRITERION_HOURS`, `ABSOLUTE_GATE_LKFS`, `RELATIVE_GATE_LU`, `SURROUND_WEIGHT` and `K_WEIGHTING_MIN_SR`) are exported in Rust and Python only; the WebAssembly binding exports functions. Errors are `SafetyError` in Rust, `ValueError` in Python (`OverflowError` for a negative or too-large integer argument, `TypeError` for a non-integer one) and `Error` in JavaScript.

### Silence, and other choices

- **Silence has no level.** A ratio of 0 is −∞ dB. Every function that returns one level rejects silence with an error: `amplitude_to_db(0)`, `rms_dbfs` of an all-zero signal, `normalize` of a signal that is silent after its high-pass, and `integrated_loudness` when every 400 ms block is below the absolute gate. For `normalize`, silent means a high-passed RMS no larger than the high-pass's rounding level, `(len(x) + L)·(ε·peak + 2⁻¹⁰⁷⁴)` for the peak of `x` and the high-pass's box length `L` (at least three times the largest residue found in an adversarial search of constants; the `2⁻¹⁰⁷⁴` term, one smallest subnormal step for each of those `len(x) + L` steps, sets the floor for subnormal samples): a constant such as `[0.3] * 1000` leaves only rounding residue of about 10⁻¹⁴, which scaling would turn into noise at the target level. `rms_dbfs` is computed relative to the peak, so tiny or huge samples (10⁻²⁰⁰ or 10³⁰⁰) have a level although their squares underflow or overflow an `f64`. `lifesound.py` added 10⁻¹² to the RMS before the logarithm, which turned silence into −240 dBFS and raised a −120 dBFS signal by 8.7 × 10⁻⁶ dB; this library does not. The two meters, which return one value per window, give −∞ for a silent window (`-inf` in Python, `-Infinity` in JavaScript), in all three languages alike.
- **K-weighting at any sample rate.** BS.1770 publishes the filter at 48 kHz only. Mapped back through the bilinear transform, its two stages are a high shelf at 1681.97 Hz (Q 0.70718, +4.00 dB above it) and a high-pass at 38.135 Hz (Q 0.50033). `k_weighting(sr)` discretises these same analogue filters at `sr`, prewarped at each corner, so at 48 kHz it returns the published coefficients again (to 10⁻¹⁴) and below 6 kHz its response at 44.1, 96 or 192 kHz stays within 0.01 dB of the 48 kHz filter. These corners and Q values are those of pyloudnorm's DeMan filters. The high-pass keeps its high-frequency gain of 1.004995 (+0.043 dB) at every rate, where pyloudnorm's DeMan filters keep its numerator `[1, −2, 1]`, so the two agree exactly at 48 kHz only. Elsewhere they differ by a nearly flat gain, which shifts loudness readings by the same amount: 0.004 dB at 44.1 kHz, 0.020 to 0.032 dB at 88.2 to 192 kHz, and more below 32 kHz (0.022 dB at 32 kHz, 0.051 at 22.05 kHz, 0.086 at 16 kHz and 0.215 at 8 kHz). The sample rate must be at least 3364 Hz, so that the shelf lies below the Nyquist frequency.
- **Blocks.** Gating blocks start every 100 ms at frame `⌊k·sr/10⌋` and are averaged over the frames they hold, so every sample rate works; at rates divisible by 10 every block is exactly 0.4·sr frames.
- **`normalize`** takes the channel count and sample rate explicitly (`normalize(x, channels, target_rms_dbfs, peak_ceiling_dbfs, sr)`), where `lifesound.py` read the channels from the array's shape and used 44.1 kHz for its 12 Hz high-pass; the high-pass is signal-core's.
- **TWA** is the exact inverse of the dose rule, `criterion + exchange · log₂(D/100)`. With a 3 dB exchange rate it is close to, but not the same as, the equal-energy Leq of the day: for 3 h at 82, 0.5 h at 97, 1.5 h at 90 and 3 h at 75 dBA it is 87.592 dBA against an Leq of 87.574 dBA. With an exchange rate of exactly 10·log₁₀ 2 dB the two agree. The form `criterion + 10·log₁₀(D/100)` treats a 3 dB exchange rate as exactly 10·log₁₀ 2 dB, so it differs from `twa_from_dose` by 0.0103 dB per doubling of the dose: at 1000 % it gives 95.0 dBA where `twa_from_dose(1000.0)` gives 94.966.
- **Zero durations, tiny ones and huge ones.** A dose segment of 0 hours counts for nothing, however loud (`noise_dose([5000, 85], [0, 8])` is 100 %). `leq` takes each segment's share of the time in dB, from the logarithms of the durations relative to the longest, so their sum may exceed the range of an `f64` and a segment 10⁻³³⁰ as long as the longest still counts (`leq([0, 6000], [1e10, 1e-320])` is about 2700 dB). `noise_dose` sums its shares `Cᵢ/Tᵢ` as powers of two around the largest, so 10⁻¹⁰ hours at 3160 dBA is a dose of 4.49 × 10²⁹⁹ %, although 2^(3075/3) overflows an `f64`; `allowed_time` then computes the time as one power of two, so 3160 dBA is allowed 2⁻¹⁰²² hours rather than 0.
- **Overflow.** Filtered signals and loudness energies are checked: `k_weight` of samples near the largest `f64`, or `integrated_loudness` of samples of about 10¹⁵² or more (whose sums of squares over a block overflow), is an error rather than ∞ or NaN. `normalize` scales the high-passed signal to a peak of 1 before its gain, so its gain overflows only when its result would not fit in an `f64`, not for a ceiling that does not bind or for subnormal samples; samples so large that the high-pass's running sums overflow, which takes a peak of at least about 1.8 × 10³⁰⁸/(len(x) + L), are an error on `x`, even when the result would fit. A channel of weight 0 takes no part in the loudness.

### Exposure limits referenced

- **NIOSH** (`safety:NIOSH-REL`, level R): the default parameters are NIOSH's recommended exposure limit, 85 dBA as an 8-hour time-weighted average with a 3 dB exchange rate (NIOSH 1998, publication 98-126, Chapter 1, Section 1.1, [doi:10.26616/NIOSHPUB98126](https://doi.org/10.26616/NIOSHPUB98126); [PDF](https://www.cdc.gov/niosh/docs/98-126/pdfs/CriteriaDoc_98-126.pdf)). The three values were checked against the document's text.
- **WHO-ITU safe listening** (WHO 2019, [Safe listening devices and systems: a WHO-ITU standard](https://www.who.int/publications/i/item/9789241515276)): reference exposures of 80 dB for 40 hours a week for adults and 75 dB for 40 hours a week for children, that is `allowed_time(L, 80.0, 3.0, 40.0)` and `allowed_time(L, 75.0, 3.0, 40.0)` in hours a week. The standard has no DOI, so it is cited here rather than as an R claim, and the 3 dB exchange rate in those calls is this library's assumption: the WHO pages that were checked do not state it, although their tables are consistent with it.
- The action values of EU Directive 2003/10/EC are not included: they could not be checked against a DOI or PMID for this release.

These are reference values for computing, not medical advice.

## How it is verified

- **Proofs:** the Lean 4 files in [`proofs/`](proofs) prove the five P-level claims in exact real arithmetic: gains add in decibels (`20·log₁₀(ab) = 20·log₁₀ a + 20·log₁₀ b`), the incoherent sum of `n` levels lies between the loudest level `M` and `M + 10·log₁₀ n`, one exchange rate more halves the allowed time (so `t` hours at `L` and `t/2` hours at `L + q` are the same dose), the TWA of a steady day returns its level, and `3.0102 < 10·log₁₀ 2 < 3.0103`. The last one reduces to the whole-number inequalities `10^30102 < 2^100000 < 10^30103`, which the Lean kernel checks. CI rebuilds the proofs with Mathlib, rejects `sorry`, and allows only the standard axioms (`propext`, `Classical.choice`, `Quot.sound`).
- **Cited, not proved:** the NIOSH limit is an R claim, checked against its document; the WHO-ITU values are cited above. BS.1770 ([ITU-R BS.1770](https://www.itu.int/rec/R-REC-BS.1770)) and EBU Tech 3341 ([tech3341.pdf](https://tech.ebu.ch/docs/tech/tech3341.pdf)) define what the loudness code is tested against; the revision of each was not checked against those pages for this release.
- **Implementation:** one Rust implementation, compiled to Python (PyO3, numpy arrays) and WebAssembly. The same 147 golden vectors in [`vectors/`](vectors) are checked in all three (109 with expected outputs from the oracle, 38 invalid inputs that must be rejected); CI runs the Rust and Python checks on Linux and Windows and the WebAssembly checks on Linux. Their oracle, [`vectors/generate.py`](vectors/generate.py), evaluates the decibel and dose formulas in 50-digit arithmetic (mpmath), runs `lifesound.py`'s `rms_dbfs` and `normalize` with numpy (with its own box filter, not signal-core's), derives the K-weighting from BS.1770's published coefficients in 50-digit arithmetic, filters with scipy, and measures loudness with pyloudnorm's meter fed the same filter, cross-checked against a direct transcription of the BS.1770 formulas. Nothing in it calls the Rust code. In Rust only, [`tests/ebu3341.rs`](tests/ebu3341.rs) measures the EBU Tech 3341 test signals 1 to 6 and 9 at 48 kHz (tens of seconds each, too long for the vector files): every reading is within the published ±0.1 LU (the meters of cases 1 and 2 within 0.007 LU, case 9 within 0.014 LU), and the integrated readings are within 10⁻⁶ LU of pyloudnorm's. Property tests in [`tests/properties.rs`](tests/properties.rs) check the proved identities in the floating-point code, the meter's scaling laws, NaN and infinite input, and finite input that overflows.
- **What the proofs don't cover:** they are about exact real numbers, not the 64-bit floating point of the code, whose rounding is checked by the golden-vector tolerances and the property tests. They prove the algebra of the dose rule, not that the rule protects hearing: the 85 dBA limit and the 3 dB exchange rate are taken from NIOSH's recommendation (a referenced value, level R, see above), not derived. Nothing about loudness is proved: K-weighting, gating and the meters are tested against BS.1770's coefficients, pyloudnorm and the EBU test signals. The EBU signals are synthesised from the test-case descriptions rather than read from EBU's files, and the programme-material and true-peak cases of Tech 3341 are not covered.

## Cite

If you use this work, please cite it. GitHub's **Cite this repository** button gives APA and BibTeX from [CITATION.cff](CITATION.cff).

## Licence

Code: Apache-2.0 ([LICENSE](LICENSE)). Copyright notices are in [NOTICE](NOTICE).
