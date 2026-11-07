// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! # Audio safety math
//!
//! Decibel arithmetic, noise-dose rules and loudness metering in Rust, Python and WebAssembly, with Lean 4 proofs of the dose rule.
//!
//! Part of [Frequencies of Life](https://chethandvg.github.io/frequencies-of-life/). Each public function corresponds to
//! claims in `claims.yaml`; the README lists what is proved and how.
//!
//! Three groups of functions:
//!
//! * **Decibels:** ratios to and from decibels ([`amplitude_to_db`] uses
//!   20·log₁₀, [`power_to_db`] 10·log₁₀), the level of unrelated sources played
//!   together ([`sum_incoherent_db`]), the equal-energy average of levels
//!   ([`energetic_mean_db`], and [`leq`] for timed segments), and the level of
//!   a signal in dB relative to full scale ([`rms_dbfs`], [`peak_dbfs`],
//!   [`normalize`], ported from the project's earlier toolkit `lifesound.py`).
//! * **Noise dose:** the exchange-rate rule ([`allowed_time`]), the daily dose
//!   ([`noise_dose`]) and the equivalent eight-hour level ([`twa_from_dose`]),
//!   with the NIOSH recommended exposure limit (85 dBA as an 8-hour TWA, 3 dB
//!   exchange rate, `safety:NIOSH-REL`) as the usual parameters.
//! * **Loudness** per ITU-R BS.1770: the K-weighting filter at any sample rate
//!   ([`k_weighting`], [`k_weight`]) and the gated [`integrated_loudness`] in
//!   LUFS, with [`momentary_loudness`] and [`short_term_loudness`].
//!
//! ```
//! use audio_safety_math::{allowed_time, integrated_loudness, sum_incoherent_db};
//! use signal_core::osc;
//!
//! // Two equal, unrelated sources are 10·log10(2) = 3.0103 dB louder than one.
//! assert!((sum_incoherent_db(&[80.0, 80.0])? - 83.0103).abs() < 1e-4);
//!
//! // NIOSH: 85 dBA is allowed 8 hours a day; every 3 dB more halves the time.
//! assert_eq!(allowed_time(88.0, 85.0, 3.0, 8.0)?, 4.0);
//! assert_eq!(allowed_time(100.0, 85.0, 3.0, 8.0)?, 0.25); // 15 minutes
//!
//! // A 1 kHz sine at -23 dBFS in both stereo channels reads -23 LUFS.
//! let tone: Vec<f64> = osc(1000.0, 2.0, 48_000, 0.0)?.iter().map(|v| 0.0707946 * v).collect();
//! let stereo: Vec<f64> = tone.iter().flat_map(|&v| [v, v]).collect();
//! let lufs = integrated_loudness(&stereo, &[1.0, 1.0], 48_000)?;
//! assert!((lufs + 23.0).abs() < 0.01);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! **Silence has no level.** A ratio of 0 is −∞ dB, so every function that
//! returns one level rejects silence (a zero ratio, an all-zero signal, a
//! signal whose every loudness block falls below the absolute gate) with a
//! [`SafetyError`] rather than return a made-up floor; `lifesound.py` added
//! 10⁻¹² before the logarithm, which turned silence into −240 dBFS. The two
//! loudness meters, which return one value per window, give −∞ for a silent
//! window. Invalid arguments (NaN, an infinite level, a negative duration, a
//! buffer that is not a whole number of frames) are a [`SafetyError`] too,
//! never a panic, and so is a result too large for an `f64`, or a filtered
//! signal or loudness energy that overflows one on the way.
//!
//! These are formulas and reference values, not medical advice or a hearing
//! test: the README says what the limits do and do not cover.

use std::fmt;

mod db;
mod dose;
mod loudness;

pub use db::{
    amplitude_to_db, db_to_amplitude, db_to_power, energetic_mean_db, leq, normalize, peak_dbfs,
    power_to_db, rms_dbfs, sum_incoherent_db,
};
pub use dose::{
    CRITERION_HOURS, NIOSH_CRITERION_DB, NIOSH_EXCHANGE_DB, allowed_time, noise_dose, twa_from_dose,
};
pub use loudness::{
    ABSOLUTE_GATE_LKFS, K_WEIGHTING_MIN_SR, RELATIVE_GATE_LU, SURROUND_WEIGHT, integrated_loudness,
    k_weight, k_weighting, momentary_loudness, short_term_loudness,
};

/// An argument (or the result) outside the range the function accepts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SafetyError {
    name: &'static str,
    value: f64,
    requirement: &'static str,
}

impl SafetyError {
    /// The parameter that was out of range (`"result"` when the result does not fit).
    pub fn name(&self) -> &'static str {
        self.name
    }

    /// The offending value (for a length or a count, the number as an `f64`).
    pub fn value(&self) -> f64 {
        self.value
    }
}

impl fmt::Display for SafetyError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} must be {}, got {}",
            self.name, self.requirement, self.value
        )
    }
}

impl std::error::Error for SafetyError {}

impl From<signal_core::SignalError> for SafetyError {
    /// Arguments are checked before signal-core sees them, so this is a fallback.
    fn from(e: signal_core::SignalError) -> Self {
        Self {
            name: e.name(),
            value: e.value(),
            requirement: "accepted by signal_core::highpass",
        }
    }
}

/// The result type of this crate.
pub type Result<T> = std::result::Result<T, SafetyError>;

const FINITE: &str = "a finite number";
const POSITIVE: &str = "a finite number greater than 0";

fn require(ok: bool, name: &'static str, value: f64, requirement: &'static str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(SafetyError {
            name,
            value,
            requirement,
        })
    }
}

fn finite(name: &'static str, x: f64) -> Result<f64> {
    require(x.is_finite(), name, x, FINITE).map(|()| x)
}

fn positive(name: &'static str, x: f64) -> Result<f64> {
    require(x.is_finite() && x > 0.0, name, x, POSITIVE).map(|()| x)
}

/// Every element of `xs` satisfies `ok`; the first that does not is the error value.
fn each(
    name: &'static str,
    xs: &[f64],
    ok: impl Fn(f64) -> bool,
    requirement: &'static str,
) -> Result<()> {
    match xs.iter().find(|&&x| !ok(x)) {
        Some(&bad) => Err(SafetyError {
            name,
            value: bad,
            requirement,
        }),
        None => Ok(()),
    }
}

fn all_finite(name: &'static str, xs: &[f64]) -> Result<()> {
    each(name, xs, f64::is_finite, "finite in every element")
}

fn all_non_negative(name: &'static str, xs: &[f64]) -> Result<()> {
    each(
        name,
        xs,
        |x| x.is_finite() && x >= 0.0,
        "finite and at least 0 in every element",
    )
}

fn same_length(name: &'static str, xs: &[f64], len: usize) -> Result<()> {
    require(
        xs.len() == len,
        name,
        xs.len() as f64,
        "as long as the list of levels",
    )
}

/// A result that must fit in an `f64`.
fn fits(x: f64) -> Result<f64> {
    require(x.is_finite(), "result", x, "small enough to fit in an f64").map(|()| x)
}

/// The number of frames in an interleaved buffer of `channels` channels.
fn frames(name: &'static str, x: &[f64], channels: usize) -> Result<usize> {
    require(channels >= 1, "channels", channels as f64, "at least 1")?;
    require(
        x.len() % channels == 0,
        name,
        x.len() as f64,
        "a whole number of frames long (a multiple of channels)",
    )?;
    Ok(x.len() / channels)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn errors_name_the_argument() {
        assert_eq!(
            frames("x", &[0.0; 3], 2).unwrap_err().to_string(),
            "x must be a whole number of frames long (a multiple of channels), got 3"
        );
        assert_eq!(frames("x", &[], 0).unwrap_err().name(), "channels");
        assert_eq!(fits(f64::INFINITY).unwrap_err().name(), "result");
        assert_eq!(
            all_finite("levels", &[1.0, f64::NAN]).unwrap_err().name(),
            "levels"
        );
        assert_eq!(
            all_non_negative("hours", &[1.0, -2.0]).unwrap_err().value(),
            -2.0
        );
    }
}
