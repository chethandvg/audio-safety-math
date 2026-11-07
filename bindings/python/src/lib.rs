// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Python bindings for `audio-safety-math`. Every function keeps the Rust name
//! and argument order, so the golden vectors can call it by name. The dose
//! functions default to the NIOSH rule (`criterion_db=85.0`, `exchange_db=3.0`,
//! `criterion_hours=8.0`); the functions that take a signal have no defaults.
//! Signals and lists of levels are accepted as numpy arrays or any sequence of
//! numbers, and signals come back as float64 numpy arrays; invalid input raises
//! `ValueError` (`OverflowError` for a negative or too-large integer argument
//! such as a channel count or sample rate, `TypeError` for a non-integer one).

use pyo3::prelude::*;

#[pymodule]
mod audio_safety_math {
    use audio_safety_math_rs as imp;
    use numpy::{AllowTypeChange, IntoPyArray, PyArray1, PyArrayLike1};
    use pyo3::exceptions::PyValueError;
    use pyo3::prelude::*;
    use std::borrow::Cow;

    /// A signal or list argument: a float64 numpy array, or anything numpy can convert to one.
    type Samples<'py> = PyArrayLike1<'py, f64, AllowTypeChange>;
    /// A signal result.
    type Array<'py> = Bound<'py, PyArray1<f64>>;

    fn ok<T>(r: imp::Result<T>) -> PyResult<T> {
        r.map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn array<'py>(py: Python<'py>, r: imp::Result<Vec<f64>>) -> PyResult<Array<'py>> {
        Ok(ok(r)?.into_pyarray(py))
    }

    /// The elements of `x`, borrowed when the array is contiguous.
    fn slice<'a>(x: &'a Samples<'_>) -> Cow<'a, [f64]> {
        x.as_slice()
            .map_or_else(|_| Cow::Owned(x.as_array().to_vec()), Cow::Borrowed)
    }

    #[pymodule_init]
    fn init(m: &Bound<'_, PyModule>) -> PyResult<()> {
        m.add("__version__", env!("CARGO_PKG_VERSION"))?;
        m.add("NIOSH_CRITERION_DB", imp::NIOSH_CRITERION_DB)?;
        m.add("NIOSH_EXCHANGE_DB", imp::NIOSH_EXCHANGE_DB)?;
        m.add("CRITERION_HOURS", imp::CRITERION_HOURS)?;
        m.add("ABSOLUTE_GATE_LKFS", imp::ABSOLUTE_GATE_LKFS)?;
        m.add("RELATIVE_GATE_LU", imp::RELATIVE_GATE_LU)?;
        m.add("SURROUND_WEIGHT", imp::SURROUND_WEIGHT)?;
        m.add("K_WEIGHTING_MIN_SR", imp::K_WEIGHTING_MIN_SR)
    }

    /// An amplitude ratio in decibels: 20 log10(ratio); ratio must be > 0.
    #[pyfunction]
    fn amplitude_to_db(ratio: f64) -> PyResult<f64> {
        ok(imp::amplitude_to_db(ratio))
    }

    /// The amplitude ratio of db decibels: 10 ** (db / 20).
    #[pyfunction]
    fn db_to_amplitude(db: f64) -> PyResult<f64> {
        ok(imp::db_to_amplitude(db))
    }

    /// A power ratio in decibels: 10 log10(ratio); ratio must be > 0.
    #[pyfunction]
    fn power_to_db(ratio: f64) -> PyResult<f64> {
        ok(imp::power_to_db(ratio))
    }

    /// The power ratio of db decibels: 10 ** (db / 10).
    #[pyfunction]
    fn db_to_power(db: f64) -> PyResult<f64> {
        ok(imp::db_to_power(db))
    }

    /// The level of unrelated sources together: 10 log10(sum(10 ** (L / 10))).
    #[pyfunction]
    fn sum_incoherent_db(levels: Samples<'_>) -> PyResult<f64> {
        ok(imp::sum_incoherent_db(&slice(&levels)))
    }

    /// The energetic (equal-energy) mean of levels: 10 log10(mean(10 ** (L / 10))).
    #[pyfunction]
    fn energetic_mean_db(levels: Samples<'_>) -> PyResult<f64> {
        ok(imp::energetic_mean_db(&slice(&levels)))
    }

    /// The equivalent continuous level of timed segments: 10 log10(sum(t 10 ** (L / 10)) / sum(t)).
    #[pyfunction]
    fn leq(levels_db: Samples<'_>, durations: Samples<'_>) -> PyResult<f64> {
        ok(imp::leq(&slice(&levels_db), &slice(&durations)))
    }

    /// The RMS level of x in dBFS: 20 log10(sqrt(mean(x ** 2))); silence raises ValueError.
    #[pyfunction]
    fn rms_dbfs(x: Samples<'_>) -> PyResult<f64> {
        ok(imp::rms_dbfs(&slice(&x)))
    }

    /// The peak level of x in dBFS: 20 log10(max(abs(x))); silence raises ValueError.
    #[pyfunction]
    fn peak_dbfs(x: Samples<'_>) -> PyResult<f64> {
        ok(imp::peak_dbfs(&slice(&x)))
    }

    /// lifesound.normalize: 12 Hz high-pass, RMS to target_rms_dbfs, peak at most peak_ceiling_dbfs.
    #[pyfunction]
    fn normalize<'py>(
        py: Python<'py>,
        x: Samples<'py>,
        channels: usize,
        target_rms_dbfs: f64,
        peak_ceiling_dbfs: f64,
        sr: u32,
    ) -> PyResult<Array<'py>> {
        array(
            py,
            imp::normalize(&slice(&x), channels, target_rms_dbfs, peak_ceiling_dbfs, sr),
        )
    }

    /// Time allowed at level_dba: criterion_hours / 2 ** ((level_dba - criterion_db) / exchange_db).
    #[pyfunction]
    #[pyo3(signature = (level_dba, criterion_db = imp::NIOSH_CRITERION_DB, exchange_db = imp::NIOSH_EXCHANGE_DB, criterion_hours = imp::CRITERION_HOURS))]
    fn allowed_time(
        level_dba: f64,
        criterion_db: f64,
        exchange_db: f64,
        criterion_hours: f64,
    ) -> PyResult<f64> {
        ok(imp::allowed_time(
            level_dba,
            criterion_db,
            exchange_db,
            criterion_hours,
        ))
    }

    /// The daily noise dose in percent: 100 * sum(hours[i] / allowed_time(levels_dba[i])).
    #[pyfunction]
    #[pyo3(signature = (levels_dba, hours, criterion_db = imp::NIOSH_CRITERION_DB, exchange_db = imp::NIOSH_EXCHANGE_DB, criterion_hours = imp::CRITERION_HOURS))]
    fn noise_dose(
        levels_dba: Samples<'_>,
        hours: Samples<'_>,
        criterion_db: f64,
        exchange_db: f64,
        criterion_hours: f64,
    ) -> PyResult<f64> {
        ok(imp::noise_dose(
            &slice(&levels_dba),
            &slice(&hours),
            criterion_db,
            exchange_db,
            criterion_hours,
        ))
    }

    /// The time-weighted average level of a dose: criterion_db + exchange_db * log2(dose_percent / 100).
    #[pyfunction]
    #[pyo3(signature = (dose_percent, criterion_db = imp::NIOSH_CRITERION_DB, exchange_db = imp::NIOSH_EXCHANGE_DB))]
    fn twa_from_dose(dose_percent: f64, criterion_db: f64, exchange_db: f64) -> PyResult<f64> {
        ok(imp::twa_from_dose(dose_percent, criterion_db, exchange_db))
    }

    /// The BS.1770 K-weighting at sr: [shelf, high_pass], each [b0, b1, b2, a1, a2] with a0 = 1.
    #[pyfunction]
    fn k_weighting(sr: u32) -> PyResult<[[f64; 5]; 2]> {
        ok(imp::k_weighting(sr))
    }

    /// x (interleaved, channels channels) through the K-weighting filter.
    #[pyfunction]
    fn k_weight<'py>(
        py: Python<'py>,
        x: Samples<'py>,
        channels: usize,
        sr: u32,
    ) -> PyResult<Array<'py>> {
        array(py, imp::k_weight(&slice(&x), channels, sr))
    }

    /// Gated integrated loudness (BS.1770) in LUFS of x, interleaved with len(weights) channels.
    #[pyfunction]
    fn integrated_loudness(x: Samples<'_>, weights: Samples<'_>, sr: u32) -> PyResult<f64> {
        ok(imp::integrated_loudness(&slice(&x), &slice(&weights), sr))
    }

    /// Momentary loudness in LUFS: every 400 ms window, one every 100 ms (-inf when silent).
    #[pyfunction]
    fn momentary_loudness<'py>(
        py: Python<'py>,
        x: Samples<'py>,
        weights: Samples<'py>,
        sr: u32,
    ) -> PyResult<Array<'py>> {
        array(
            py,
            imp::momentary_loudness(&slice(&x), &slice(&weights), sr),
        )
    }

    /// Short-term loudness in LUFS: every 3 s window, one every 100 ms (-inf when silent).
    #[pyfunction]
    fn short_term_loudness<'py>(
        py: Python<'py>,
        x: Samples<'py>,
        weights: Samples<'py>,
        sr: u32,
    ) -> PyResult<Array<'py>> {
        array(
            py,
            imp::short_term_loudness(&slice(&x), &slice(&weights), sr),
        )
    }
}
