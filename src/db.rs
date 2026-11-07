// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Decibel arithmetic: ratios, sums and averages of levels, and signal levels in dBFS.

use crate::{
    Result, all_finite, all_non_negative, each, finite, fits, frames, positive, require,
    same_length,
};

const RATIO: &str =
    "a finite number greater than 0 (a ratio of 0, silence, has no level in decibels)";
const SIGNAL: &str =
    "a signal that is not silent: a finite level above 0 (silence has no level in decibels)";
const HIGH_PASSED: &str = "a signal that is not silent after the 12 Hz high-pass: its RMS there, in dBFS, above the high-pass's rounding level 20·log10((len(x) + L)·(ε·peak + 2^-1074)), L the box length (a constant is silent)";
const LEVELS: &str = "at least one level";

/// The cutoff of the high-pass [`normalize`] applies first, as in `lifesound.py`.
const NORMALIZE_HIGHPASS_HZ: f64 = 12.0;

/// `scale · log₁₀ r` for a ratio `r > 0`: 20 for amplitudes, 10 for powers.
fn to_db(name: &'static str, r: f64, scale: f64, requirement: &'static str) -> Result<f64> {
    require(r.is_finite() && r > 0.0, name, r, requirement)?;
    Ok(scale * r.log10())
}

/// `10^(db / scale)`; it may round to 0, but must not overflow.
fn from_db(db: f64, scale: f64) -> Result<f64> {
    fits(10f64.powf(finite("db", db)? / scale))
}

/// An amplitude (field) ratio in decibels: `20 · log₁₀ ratio`. A gain of 2 is
/// +6.02 dB; gains in a chain multiply, so their decibels add (proved as
/// `safety:DB-PRODUCT`).
pub fn amplitude_to_db(ratio: f64) -> Result<f64> {
    to_db("ratio", ratio, 20.0, RATIO)
}

/// The amplitude ratio of `db` decibels: `10^(db / 20)`.
pub fn db_to_amplitude(db: f64) -> Result<f64> {
    from_db(db, 20.0)
}

/// A power (energy, intensity) ratio in decibels: `10 · log₁₀ ratio`. Doubling
/// the power adds 10·log₁₀ 2 ≈ 3.0103 dB (`safety:THREE-DB-IS-ABOUT-DOUBLING`).
pub fn power_to_db(ratio: f64) -> Result<f64> {
    to_db("ratio", ratio, 10.0, RATIO)
}

/// The power ratio of `db` decibels: `10^(db / 10)`.
pub fn db_to_power(db: f64) -> Result<f64> {
    from_db(db, 10.0)
}

/// `10 · log₁₀ Σ 10^((Lᵢ + gᵢ)/10)` over the levels whose gain `gᵢ` (in dB)
/// is above −∞, factored around the largest term, level `M` with gain `G`, as
/// `M + (G + 10 · log₁₀ Σ 10^(((Lᵢ − M) + (gᵢ − G))/10))`, so that no term
/// overflows or underflows on the way, and gains of 0 add no rounding.
fn power_sum_db(levels: &[f64], gain_db: impl Fn(usize) -> f64) -> Result<f64> {
    let terms = || {
        (0..levels.len())
            .map(|i| (levels[i], gain_db(i)))
            .filter(|&(_, g)| g > f64::NEG_INFINITY)
    };
    let (top, top_gain) = terms().fold((f64::NEG_INFINITY, 0.0), |m, t| {
        if t.0 + t.1 > m.0 + m.1 { t } else { m }
    });
    let sum: f64 = terms()
        .map(|(l, g)| 10f64.powf(((l - top) + (g - top_gain)) / 10.0))
        .sum();
    fits(top + (top_gain + 10.0 * sum.log10()))
}

fn checked_levels(name: &'static str, levels: &[f64]) -> Result<f64> {
    require(!levels.is_empty(), name, 0.0, LEVELS)?;
    all_finite(name, levels)?;
    Ok(levels.len() as f64)
}

/// The level of unrelated (incoherent) sources played together:
/// `10 · log₁₀ Σ 10^(Lᵢ/10)`, since their powers add. Two equal sources are
/// 3.01 dB louder than one; the sum lies between the loudest level and that
/// level plus `10 · log₁₀ n` (proved as `safety:INCOHERENT-SUM-BOUNDS`).
pub fn sum_incoherent_db(levels: &[f64]) -> Result<f64> {
    checked_levels("levels", levels)?;
    power_sum_db(levels, |_| 0.0)
}

/// The energetic (equal-energy) mean of levels: `10 · log₁₀ ((1/n) Σ 10^(Lᵢ/10))`,
/// the level of the average power, which is [`sum_incoherent_db`] minus
/// `10 · log₁₀ n`. It is never below the arithmetic mean of the levels.
pub fn energetic_mean_db(levels: &[f64]) -> Result<f64> {
    let n = checked_levels("levels", levels)?;
    let share = -10.0 * n.log10();
    power_sum_db(levels, |_| share)
}

/// The equivalent continuous level of timed segments (equal energy):
/// `Leq = 10 · log₁₀ (Σ tᵢ · 10^(Lᵢ/10) / Σ tᵢ)` for segment `i` at level
/// `levels_db[i]` for `durations[i]` (any one unit). Segments of zero duration
/// count for nothing; at least one must last longer than 0. Each segment's
/// share of the time, `tᵢ / Σ tⱼ`, enters in dB, from the logarithms of the
/// durations relative to the longest, so their sum may exceed the range of an
/// `f64`, and a segment far shorter than the longest still counts.
pub fn leq(levels_db: &[f64], durations: &[f64]) -> Result<f64> {
    checked_levels("levels_db", levels_db)?;
    same_length("durations", durations, levels_db.len())?;
    all_non_negative("durations", durations)?;
    let longest = positive("durations", durations.iter().copied().fold(0.0, f64::max))?;
    let total_db = 10.0 * durations.iter().map(|t| t / longest).sum::<f64>().log10();
    let longest_log = longest.log10();
    // log10(0) = −∞, so a segment of zero duration has a gain of −∞ and takes no part.
    power_sum_db(levels_db, |i| {
        10.0 * (durations[i].log10() - longest_log) - total_db
    })
}

fn checked_signal(x: &[f64]) -> Result<()> {
    require(!x.is_empty(), "x", 0.0, "at least one sample")?;
    all_finite("x", x)
}

/// The RMS level of `x` in dB relative to full scale (an amplitude of 1):
/// `20 · log₁₀ √(mean(x²))`, over every sample of every channel. It is computed
/// relative to the peak, so the mean square neither underflows nor overflows
/// for any finite samples. `lifesound.py` added 10⁻¹² to the RMS before the
/// logarithm; this does not, so silence is an error instead of −240 dBFS.
pub fn rms_dbfs(x: &[f64]) -> Result<f64> {
    checked_signal(x)?;
    let level = rms_db(x);
    require(level.is_finite(), "x", 0.0, SIGNAL).map(|()| level)
}

/// The peak level of `x` in dB relative to full scale: `20 · log₁₀ max |xᵢ|`.
pub fn peak_dbfs(x: &[f64]) -> Result<f64> {
    checked_signal(x)?;
    to_db("x", peak(x), 20.0, SIGNAL)
}

fn peak(x: &[f64]) -> f64 {
    x.iter().fold(0.0, |m, v| m.max(v.abs()))
}

/// The peak `p` of `x` and the RMS level relative to it,
/// `10 · log₁₀ mean((x/p)²)`, whose mean square lies in `[1/n, 1]`, so it
/// neither underflows nor overflows; `(0, −∞)` for silence.
fn peak_and_relative_rms(x: &[f64]) -> (f64, f64) {
    let top = peak(x);
    if top == 0.0 {
        return (0.0, f64::NEG_INFINITY);
    }
    let mean_square = x.iter().map(|v| (v / top).powi(2)).sum::<f64>() / x.len() as f64;
    (top, 10.0 * mean_square.log10())
}

/// `20 · log₁₀ √(mean(x²))` as `20 · log₁₀ p + 10 · log₁₀ mean((x/p)²)` for the
/// peak `p`; −∞ for silence.
fn rms_db(x: &[f64]) -> f64 {
    let (top, relative) = peak_and_relative_rms(x);
    20.0 * top.log10() + relative
}

/// `lifesound.normalize`: remove DC and sub-audio drift with signal-core's
/// zero-phase 12 Hz [`highpass`](signal_core::highpass), scale to an RMS of
/// `target_rms_dbfs`, or, if a sample would then exceed `peak_ceiling_dbfs`,
/// scale so that the peak sits at the ceiling (the RMS is then lower than the
/// target). `x` is interleaved with `channels` channels and the RMS is taken
/// over all of them.
///
/// A signal that is silent after the high-pass is an error: all zeros, a
/// constant, or anything whose high-passed RMS is no more than the rounding
/// level of the high-pass's running sums, `(len(x) + L) · (ε · peak + 2⁻¹⁰⁷⁴)`
/// for the peak of `x` and the box length `L` (at least three times the
/// largest residue found in an adversarial search of constants; for subnormal
/// samples, at least `len(x) + L` of the smallest subnormal steps). Scaling
/// that residue would only turn rounding error into noise at the target level.
/// The high-passed signal is scaled to a peak of 1 before the gain is applied,
/// so the gain overflows only when the result would not fit in an `f64`,
/// however high a ceiling that does not bind, or however small (even subnormal)
/// the samples. Samples so large that the high-pass's running sums overflow,
/// which takes a peak of at least about 1.8 × 10³⁰⁸ / (len(x) + L), are an
/// error on `x`, even when the result would fit.
pub fn normalize(
    x: &[f64],
    channels: usize,
    target_rms_dbfs: f64,
    peak_ceiling_dbfs: f64,
    sr: u32,
) -> Result<Vec<f64>> {
    frames("x", x, channels)?;
    checked_signal(x)?;
    let target = finite("target_rms_dbfs", target_rms_dbfs)?;
    let ceiling = finite("peak_ceiling_dbfs", peak_ceiling_dbfs)?;
    require(sr > 0, "sr", 0.0, "at least 1 sample per second")?;
    let mut y = signal_core::highpass(x, channels, NORMALIZE_HIGHPASS_HZ, sr)?;
    let fit = "small enough that its 12 Hz high-pass fits in an f64";
    each("x", &y, f64::is_finite, fit)?;
    // The rounding level of the high-pass's running sums: (len(x) + L)·ε times
    // the peak of x, and never less than that many of the smallest subnormal steps.
    let steps = (x.len() + signal_core::box_length(NORMALIZE_HIGHPASS_HZ, sr)?) as f64;
    let rounding = steps * (f64::EPSILON * peak(x) + f64::from_bits(1));
    let (top, relative) = peak_and_relative_rms(&y);
    let rms = 20.0 * top.log10() + relative;
    require(rms > 20.0 * rounding.log10(), "x", rms, HIGH_PASSED)?;
    // Scaled to a peak of 1, the signal's RMS is `relative` dBFS; the gain in dB
    // then reaches the target, or puts the peak at the ceiling, and is the
    // output's peak level, so it overflows only when the output would.
    let gain = db_to_amplitude((target - relative).min(ceiling))?;
    y.iter_mut().for_each(|v| *v = *v / top * gain);
    Ok(y)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ratios_round_trip() {
        assert_eq!(amplitude_to_db(10.0).unwrap(), 20.0);
        assert_eq!(power_to_db(100.0).unwrap(), 20.0);
        assert_eq!(db_to_amplitude(-20.0).unwrap(), 0.1);
        assert_eq!(db_to_power(30.0).unwrap(), 1000.0);
        assert_eq!(db_to_power(-4000.0).unwrap(), 0.0); // rounds to 0
    }

    #[test]
    fn silence_and_overflow_are_errors() {
        assert_eq!(amplitude_to_db(0.0).unwrap_err().name(), "ratio");
        assert!(power_to_db(-1.0).is_err());
        assert!(power_to_db(f64::INFINITY).is_err());
        assert_eq!(db_to_power(4000.0).unwrap_err().name(), "result");
        assert_eq!(db_to_amplitude(f64::NAN).unwrap_err().name(), "db");
        assert_eq!(
            rms_dbfs(&[0.0; 4]).unwrap_err().to_string(),
            "x must be a signal that is not silent: a finite level above 0 (silence has no level in decibels), got 0"
        );
        assert!(rms_dbfs(&[]).is_err());
        assert!(peak_dbfs(&[0.0]).is_err());
        assert!(normalize(&[0.25; 64], 1, -20.0, -1.0, 1000).is_err()); // a constant is silent after the high-pass
    }

    #[test]
    fn signal_levels_have_no_underflow_or_overflow() {
        // x² would underflow to 0 or overflow to ∞; the level relative to the peak does neither.
        assert!((rms_dbfs(&[1e-200; 4]).unwrap() + 4000.0).abs() < 1e-9);
        assert!((rms_dbfs(&[1e300, -1e300]).unwrap() - 6000.0).abs() < 1e-9);
        assert!((rms_dbfs(&[5e-324]).unwrap() - peak_dbfs(&[5e-324]).unwrap()).abs() < 1e-9);
        // A high-pass that overflows is an error, not NaN.
        let huge: Vec<f64> = (0..100).map(|i| 1e307 * (1.5 + f64::from(i % 2))).collect();
        assert_eq!(
            normalize(&huge, 1, -20.0, -1.0, 1000).unwrap_err().name(),
            "x"
        );
    }

    /// Constants of every size, length and rate leave only rounding residue
    /// after the high-pass, and every one is rejected rather than scaled up.
    #[test]
    fn constants_are_silent_after_the_high_pass() {
        let mut rng = signal_core::SplitMix64::new(9);
        for n in [1, 2, 3, 10, 64, 999, 1000, 4410, 20_000] {
            for sr in [1, 50, 1000, 8000, 44_100, 48_000, 192_000] {
                for channels in [1, 2] {
                    let c: Vec<f64> = (0..channels)
                        .map(|_| rng.uniform(-1.0, 1.0) * 10f64.powf(rng.uniform(-300.0, 300.0)))
                        .collect();
                    let x: Vec<f64> = c.iter().copied().cycle().take(n * channels).collect();
                    let err = normalize(&x, channels, -20.0, -1.0, sr).unwrap_err();
                    assert_eq!(err.name(), "x", "{n} {sr} {c:?}");
                }
            }
        }
        for (c, sr) in [(0.3, 1000), (0.1, 8000), (0.123456, 48_000), (1e-310, 1000)] {
            assert!(normalize(&[c; 20_000], 1, -20.0, -1.0, sr).is_err());
        }
        // A tone on a large offset is not silent.
        let tone: Vec<f64> = (0..2000)
            .map(|i| 0.3 + 0.01 * (0.2 * f64::from(i)).sin())
            .collect();
        let y = normalize(&tone, 1, -20.0, -1.0, 1000).unwrap();
        assert!((rms_dbfs(&y).unwrap() + 20.0).abs() < 1e-9);
    }

    /// The gain overflows only for a result too large for an f64: not for a
    /// ceiling that does not bind, nor to lift subnormal samples.
    #[test]
    fn normalize_gain_overflows_only_when_the_result_does() {
        let tone = |amp: f64| -> Vec<f64> {
            (0..1000)
                .map(|i| amp * (0.2 * f64::from(i)).sin())
                .collect()
        };
        for ceiling in [0.0, 6000.0, 7000.0, 1e300] {
            let y = normalize(&tone(0.5), 1, -20.0, ceiling, 1000).unwrap();
            assert!((rms_dbfs(&y).unwrap() + 20.0).abs() < 1e-9, "{ceiling}");
        }
        for amp in [1e-300, 1e-308, 1e-310] {
            let y = normalize(&tone(amp), 1, -20.0, -1.0, 1000).unwrap();
            assert!((rms_dbfs(&y).unwrap() + 20.0).abs() < 1e-9, "{amp}");
        }
        // The ceiling binds: the peak sits at it, here 6000 dBFS (10^300).
        let y = normalize(&tone(0.5), 1, 7000.0, 6000.0, 1000).unwrap();
        assert!((peak_dbfs(&y).unwrap() - 6000.0).abs() < 1e-9);
        // A target RMS of 7000 dBFS does not fit in an f64.
        let err = normalize(&tone(0.5), 1, 7000.0, 7000.0, 1000).unwrap_err();
        assert_eq!(err.name(), "result");
    }

    #[test]
    fn sums_stay_finite_for_extreme_levels() {
        assert_eq!(
            sum_incoherent_db(&[400.0, 400.0]).unwrap(),
            400.0 + 10.0 * 2f64.log10()
        );
        assert_eq!(sum_incoherent_db(&[-5000.0]).unwrap(), -5000.0);
        assert_eq!(
            energetic_mean_db(&[7000.0, 7000.0, 7000.0]).unwrap(),
            7000.0
        );
        // A zero-length segment does not take part, however loud it is.
        assert_eq!(leq(&[1e300, 60.0], &[0.0, 2.0]).unwrap(), 60.0);
        // Durations whose sum overflows an f64 are still valid.
        assert_eq!(leq(&[60.0, 60.0], &[1e308, 1e308]).unwrap(), 60.0);
        // A segment 10^-330 as long as the longest still counts: 6000 - 3300 dB
        // (to 1e-4, since the subnormal nearest 1e-320 is about 1.1e-5 short of it).
        let tiny = leq(&[0.0, 6000.0], &[1e10, 1e-320]).unwrap();
        assert!((tiny - 2700.0).abs() < 1e-4, "{tiny}");
        assert!(sum_incoherent_db(&[]).is_err());
        assert!(leq(&[60.0], &[0.0]).is_err());
        assert!(leq(&[60.0, 70.0], &[1.0]).is_err());
        assert!(leq(&[60.0], &[-1.0]).is_err());
    }
}
