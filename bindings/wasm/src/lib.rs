// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! WebAssembly bindings for `audio-safety-math`. Every function keeps the Rust name, so
//! the golden vectors can call it by name; invalid input throws a JS `Error`.
//! Signals and lists of levels cross the boundary as `Float64Array` (plain
//! arrays of numbers are accepted too). Channel counts and sample rates are
//! ordinary JS numbers that must be whole, so JavaScript never sees a `BigInt`.
//! The dose functions take the NIOSH rule (85 dB, 3 dB, 8 h) for any trailing
//! argument left `undefined`.

use audio_safety_math_rs as imp;
use wasm_bindgen::prelude::*;

/// The largest whole number a JS number holds exactly, 2⁵³ − 1.
const MAX_SAFE_INTEGER: f64 = 9_007_199_254_740_991.0;

/// `x` as a whole number of type `T`, or an error naming the parameter.
fn whole<T: TryFrom<u64>>(name: &str, x: f64) -> Result<T, JsError> {
    if (0.0..=MAX_SAFE_INTEGER).contains(&x)
        && x.fract() == 0.0
        && let Ok(v) = T::try_from(x as u64)
    {
        return Ok(v);
    }
    Err(JsError::new(&format!(
        "{name} must be a whole number in range, got {x}"
    )))
}

/// The crate version.
#[wasm_bindgen]
pub fn version() -> String {
    env!("CARGO_PKG_VERSION").to_owned()
}

/// An amplitude ratio in decibels: 20·log₁₀ ratio; `ratio` must be > 0.
#[wasm_bindgen]
pub fn amplitude_to_db(ratio: f64) -> Result<f64, JsError> {
    Ok(imp::amplitude_to_db(ratio)?)
}

/// The amplitude ratio of `db` decibels: 10^(db / 20).
#[wasm_bindgen]
pub fn db_to_amplitude(db: f64) -> Result<f64, JsError> {
    Ok(imp::db_to_amplitude(db)?)
}

/// A power ratio in decibels: 10·log₁₀ ratio; `ratio` must be > 0.
#[wasm_bindgen]
pub fn power_to_db(ratio: f64) -> Result<f64, JsError> {
    Ok(imp::power_to_db(ratio)?)
}

/// The power ratio of `db` decibels: 10^(db / 10).
#[wasm_bindgen]
pub fn db_to_power(db: f64) -> Result<f64, JsError> {
    Ok(imp::db_to_power(db)?)
}

/// The level of unrelated sources together: 10·log₁₀ Σ 10^(L/10).
#[wasm_bindgen]
pub fn sum_incoherent_db(levels: &[f64]) -> Result<f64, JsError> {
    Ok(imp::sum_incoherent_db(levels)?)
}

/// The energetic (equal-energy) mean of levels.
#[wasm_bindgen]
pub fn energetic_mean_db(levels: &[f64]) -> Result<f64, JsError> {
    Ok(imp::energetic_mean_db(levels)?)
}

/// The equivalent continuous level of timed segments: 10·log₁₀(Σ tᵢ·10^(Lᵢ/10) / Σ tᵢ).
#[wasm_bindgen]
pub fn leq(levels_db: &[f64], durations: &[f64]) -> Result<f64, JsError> {
    Ok(imp::leq(levels_db, durations)?)
}

/// The RMS level of `x` in dBFS; silence throws.
#[wasm_bindgen]
pub fn rms_dbfs(x: &[f64]) -> Result<f64, JsError> {
    Ok(imp::rms_dbfs(x)?)
}

/// The peak level of `x` in dBFS; silence throws.
#[wasm_bindgen]
pub fn peak_dbfs(x: &[f64]) -> Result<f64, JsError> {
    Ok(imp::peak_dbfs(x)?)
}

/// `lifesound.normalize`: 12 Hz high-pass, RMS to `target_rms_dbfs`, peak at most `peak_ceiling_dbfs`.
#[wasm_bindgen]
pub fn normalize(
    x: &[f64],
    channels: f64,
    target_rms_dbfs: f64,
    peak_ceiling_dbfs: f64,
    sr: f64,
) -> Result<Vec<f64>, JsError> {
    let (channels, sr) = (whole("channels", channels)?, whole("sr", sr)?);
    Ok(imp::normalize(
        x,
        channels,
        target_rms_dbfs,
        peak_ceiling_dbfs,
        sr,
    )?)
}

/// Time allowed at `level_dba`: criterion_hours / 2^((level_dba − criterion_db) / exchange_db).
#[wasm_bindgen]
pub fn allowed_time(
    level_dba: f64,
    criterion_db: Option<f64>,
    exchange_db: Option<f64>,
    criterion_hours: Option<f64>,
) -> Result<f64, JsError> {
    Ok(imp::allowed_time(
        level_dba,
        criterion_db.unwrap_or(imp::NIOSH_CRITERION_DB),
        exchange_db.unwrap_or(imp::NIOSH_EXCHANGE_DB),
        criterion_hours.unwrap_or(imp::CRITERION_HOURS),
    )?)
}

/// The daily noise dose in percent: 100·Σ hours[i] / allowed_time(levels_dba[i]).
#[wasm_bindgen]
pub fn noise_dose(
    levels_dba: &[f64],
    hours: &[f64],
    criterion_db: Option<f64>,
    exchange_db: Option<f64>,
    criterion_hours: Option<f64>,
) -> Result<f64, JsError> {
    Ok(imp::noise_dose(
        levels_dba,
        hours,
        criterion_db.unwrap_or(imp::NIOSH_CRITERION_DB),
        exchange_db.unwrap_or(imp::NIOSH_EXCHANGE_DB),
        criterion_hours.unwrap_or(imp::CRITERION_HOURS),
    )?)
}

/// The time-weighted average level of a dose: criterion_db + exchange_db·log₂(dose / 100).
#[wasm_bindgen]
pub fn twa_from_dose(
    dose_percent: f64,
    criterion_db: Option<f64>,
    exchange_db: Option<f64>,
) -> Result<f64, JsError> {
    Ok(imp::twa_from_dose(
        dose_percent,
        criterion_db.unwrap_or(imp::NIOSH_CRITERION_DB),
        exchange_db.unwrap_or(imp::NIOSH_EXCHANGE_DB),
    )?)
}

/// The BS.1770 K-weighting at `sr`: `[shelf, highPass]`, each a `Float64Array` `[b0, b1, b2, a1, a2]`.
#[wasm_bindgen]
pub fn k_weighting(sr: f64) -> Result<Vec<JsValue>, JsError> {
    let stages = imp::k_weighting(whole("sr", sr)?)?;
    Ok(stages.iter().map(|s| s.to_vec().into()).collect())
}

/// `x` (interleaved, `channels` channels) through the K-weighting filter.
#[wasm_bindgen]
pub fn k_weight(x: &[f64], channels: f64, sr: f64) -> Result<Vec<f64>, JsError> {
    let (channels, sr) = (whole("channels", channels)?, whole("sr", sr)?);
    Ok(imp::k_weight(x, channels, sr)?)
}

/// Gated integrated loudness (BS.1770) in LUFS of `x`, interleaved with `weights.length` channels.
#[wasm_bindgen]
pub fn integrated_loudness(x: &[f64], weights: &[f64], sr: f64) -> Result<f64, JsError> {
    Ok(imp::integrated_loudness(x, weights, whole("sr", sr)?)?)
}

/// Momentary loudness in LUFS: every 400 ms window, one every 100 ms (−Infinity when silent).
#[wasm_bindgen]
pub fn momentary_loudness(x: &[f64], weights: &[f64], sr: f64) -> Result<Vec<f64>, JsError> {
    Ok(imp::momentary_loudness(x, weights, whole("sr", sr)?)?)
}

/// Short-term loudness in LUFS: every 3 s window, one every 100 ms (−Infinity when silent).
#[wasm_bindgen]
pub fn short_term_loudness(x: &[f64], weights: &[f64], sr: f64) -> Result<Vec<f64>, JsError> {
    Ok(imp::short_term_loudness(x, weights, whole("sr", sr)?)?)
}
