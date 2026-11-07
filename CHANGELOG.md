# Changelog

All notable changes are listed here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-11-07

### Added

- Decibel arithmetic: `amplitude_to_db`, `db_to_amplitude`, `power_to_db`, `db_to_power`, `sum_incoherent_db` (factored around the loudest level, so extreme levels stay finite), `energetic_mean_db` and `leq` for timed segments (the durations enter as logarithms, so a huge total or a tiny segment stays in range).
- Signal levels ported from the earlier `lifesound.py` toolkit: `rms_dbfs` (computed relative to the peak, so it neither underflows nor overflows), `peak_dbfs` and `normalize` (on signal-core's zero-phase 12 Hz `highpass`), without `lifesound.py`'s 10⁻¹² before the logarithm, so silence is an error rather than −240 dBFS; `normalize` also rejects a signal left with only rounding residue after its high-pass, such as a constant, and scales to a peak of 1 before its gain, so its gain overflows only when its result would (samples large enough to overflow the high-pass are an error).
- Noise dose by the exchange-rate rule: `allowed_time` (one power of two when `2^((L − L₀)/q)` leaves the normal `f64` range, so a time that fits is not lost), `noise_dose` (percent, its shares summed as powers of two around the largest, so a short time at a huge level stays in range) and `twa_from_dose`, with the NIOSH recommended exposure limit (level R, checked against NIOSH 1998, Section 1.1) in `NIOSH_CRITERION_DB` (85 dBA), `NIOSH_EXCHANGE_DB` (3 dB) and `CRITERION_HOURS` (8 h).
- Loudness per ITU-R BS.1770: `k_weighting` at any sample rate from 3364 Hz (the published 48 kHz filter mapped to analogue prototypes and back by the prewarped bilinear transform), `k_weight`, the gated `integrated_loudness` in LUFS, and `momentary_loudness` and `short_term_loudness`; constants `ABSOLUTE_GATE_LKFS`, `RELATIVE_GATE_LU`, `SURROUND_WEIGHT` and `K_WEIGHTING_MIN_SR`.
- `SafetyError` for invalid input and for results that do not fit in an `f64`, including filtered signals and loudness energies that overflow.
- Python (PyO3, numpy arrays, NIOSH defaults for the dose functions) and WebAssembly (`Float64Array`, NIOSH values for omitted dose arguments) bindings with the same function names.
- Lean 4 proofs of five P-level claims: decibels of products add, the bounds of an incoherent sum, the exchange-rate halving and its equal doses, the TWA of a steady day, and 3.0102 < 10·log₁₀ 2 < 3.0103; the NIOSH limit is an R claim, checked against NIOSH 1998, Section 1.1 (doi:10.26616/NIOSHPUB98126).
- Golden vectors (147 cases) shared by Rust, Python and WebAssembly: 109 whose expected outputs come from mpmath, `lifesound.py` with numpy, scipy, pyloudnorm or a transcription of BS.1770, and 38 invalid inputs, from the crate's own input rules, that must be rejected; the EBU Tech 3341 test signals 1 to 6 and 9 at 48 kHz (in Rust); and property tests.
