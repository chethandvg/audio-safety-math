// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! EBU Tech 3341 conformance at 48 kHz: the test signals described in
//! `tests/ebu3341.json` (written by `vectors/generate.py`), measured here and
//! compared with the published expected values (±0.1 LU) and with pyloudnorm's
//! integrated loudness of the same signals. The signals are tens of seconds
//! long, too long for the shared vector files, so this check runs in Rust only;
//! the bindings call the same code.

use audio_safety_math::{
    db_to_amplitude, integrated_loudness, momentary_loudness, short_term_loudness,
};
use serde_json::Value;

/// Agreement with pyloudnorm: its DeMan filters reproduce the published
/// coefficients to 10⁻¹², and both meters sum the same blocks in double
/// precision, so 10⁻⁶ LU leaves a wide margin for rounding over millions of samples.
const PYLOUDNORM_LU: f64 = 1e-6;

/// The tighter bounds that `safety:EBU-3341` states for the meters, inside the
/// published ±0.1 LU: cases 1 and 2 within 0.007 LU of the expected value,
/// case 9 within 0.014 LU.
fn claimed_tolerance(case: u64, measure: &str, published: f64) -> f64 {
    match (case, measure) {
        (1 | 2, "momentary" | "short_term") => 0.007,
        (9, "short_term") => 0.014,
        _ => published,
    }
}

fn floats(v: &Value) -> Vec<f64> {
    v.as_array()
        .expect("an array")
        .iter()
        .map(|x| x.as_f64().expect("a number"))
        .collect()
}

/// The case's 1 kHz sine, interleaved, each segment at its level in dBFS
/// (one level for every channel, or one per channel).
fn signal(sr: u32, hz: f64, channels: usize, segments: &[Value]) -> Vec<f64> {
    let parts: Vec<(usize, Vec<f64>)> = segments
        .iter()
        .map(|s| {
            let n = (s[0].as_f64().unwrap() * f64::from(sr)).round() as usize;
            let level = &s[1];
            let dbfs = if level.is_array() {
                floats(level)
            } else {
                vec![level.as_f64().unwrap(); channels]
            };
            (
                n,
                dbfs.iter().map(|&d| db_to_amplitude(d).unwrap()).collect(),
            )
        })
        .collect();
    let total: usize = parts.iter().map(|(n, _)| n).sum();
    let tone = signal_core::osc(hz, total as f64 / f64::from(sr), sr, 0.0).unwrap();
    assert_eq!(tone.len(), total);
    let gains = parts.iter().flat_map(|(n, g)| std::iter::repeat_n(g, *n));
    tone.iter()
        .zip(gains)
        .flat_map(|(&v, g)| g.iter().map(move |a| a * v))
        .collect()
}

#[test]
fn ebu_tech_3341_test_signals() {
    let doc: Value = serde_json::from_str(include_str!("ebu3341.json")).unwrap();
    let sr = u32::try_from(doc["sr"].as_u64().unwrap()).unwrap();
    let hz = doc["hz"].as_f64().unwrap();
    for case in doc["cases"].as_array().unwrap() {
        let weights = floats(&case["weights"]);
        let x = signal(sr, hz, weights.len(), case["segments"].as_array().unwrap());
        let (expected, tol) = (
            case["expected"].as_f64().unwrap(),
            case["tolerance"].as_f64().unwrap(),
        );
        let (number, measure) = (&case["case"], case["measure"].as_str().unwrap());
        let tol = claimed_tolerance(number.as_u64().unwrap(), measure, tol);
        let values = match measure {
            "integrated" => {
                let got = integrated_loudness(&x, &weights, sr).unwrap();
                let reference = case["pyloudnorm"].as_f64().unwrap();
                assert!(
                    (got - reference).abs() < PYLOUDNORM_LU,
                    "case {number}: {got} LUFS, pyloudnorm {reference}"
                );
                vec![got]
            }
            "momentary" => momentary_loudness(&x, &weights, sr).unwrap(),
            "short_term" => short_term_loudness(&x, &weights, sr).unwrap(),
            other => panic!("unknown measure {other}"),
        };
        let (lo, hi) = values
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
                (lo.min(v), hi.max(v))
            });
        println!(
            "EBU Tech 3341 case {number}, {measure}: {lo:.4} to {hi:.4} LUFS over {} values (expected {expected} ± {tol})",
            values.len()
        );
        assert!(
            !values.is_empty() && expected - tol <= lo && hi <= expected + tol,
            "case {number} {measure}: {lo} to {hi}, expected {expected} ± {tol}"
        );
    }
}
