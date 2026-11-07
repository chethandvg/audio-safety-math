// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Property tests over many inputs, for what single golden vectors cannot show:
//! the proved identities in the floating-point code, scaling laws of the
//! loudness meter, and NaN or infinite input (which JSON cannot hold).

use audio_safety_math::*;
use signal_core::{SplitMix64, osc};
use std::f64::consts::PI;

/// `n` pseudo-random levels in [lo, hi) dB.
fn levels(rng: &mut SplitMix64, n: usize, lo: f64, hi: f64) -> Vec<f64> {
    (0..n).map(|_| rng.uniform(lo, hi)).collect()
}

/// safety:INCOHERENT-SUM-BOUNDS in the code: max ≤ sum ≤ max + 10·log10 n, exactly, in f64.
#[test]
fn incoherent_sum_stays_within_the_proved_bounds() {
    let mut rng = SplitMix64::new(1);
    for n in 1..=40 {
        for (lo, hi) in [(-50.0, 120.0), (59.0, 61.0), (-4000.0, 4000.0)] {
            let l = levels(&mut rng, n, lo, hi);
            let top = l.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let sum = sum_incoherent_db(&l).unwrap();
            assert!(
                top <= sum && sum <= top + 10.0 * (n as f64).log10(),
                "{l:?}"
            );
            let mean = energetic_mean_db(&l).unwrap();
            let arithmetic = l.iter().sum::<f64>() / n as f64;
            assert!((mean - (sum - 10.0 * (n as f64).log10())).abs() < 1e-9);
            assert!(mean >= arithmetic - 1e-9 && mean <= top);
            let equal = leq(&l, &vec![2.5; n]).unwrap();
            assert!((equal - mean).abs() <= 1e-12 * mean.abs().max(1.0));
        }
    }
}

/// safety:DB-PRODUCT and safety:THREE-DB-IS-ABOUT-DOUBLING in the code.
#[test]
fn decibels_of_products_add_and_doubling_is_three_db() {
    let mut rng = SplitMix64::new(2);
    for _ in 0..2000 {
        let (a, b) = (
            10f64.powf(rng.uniform(-6.0, 6.0)),
            10f64.powf(rng.uniform(-6.0, 6.0)),
        );
        let sum = amplitude_to_db(a).unwrap() + amplitude_to_db(b).unwrap();
        assert!((amplitude_to_db(a * b).unwrap() - sum).abs() < 1e-12);
        let back = db_to_amplitude(amplitude_to_db(a).unwrap()).unwrap();
        assert!((back - a).abs() <= 1e-13 * a);
    }
    let doubling = power_to_db(2.0).unwrap();
    assert!(3.0102 < doubling && doubling < 3.0103);
    assert_eq!(sum_incoherent_db(&[70.0, 70.0]).unwrap(), 70.0 + doubling);
}

/// safety:EXCHANGE-RATE and safety:TWA-INVERTS-DOSE in the code.
#[test]
fn one_exchange_rate_halves_the_time_and_twa_inverts_the_dose() {
    let mut rng = SplitMix64::new(3);
    for _ in 0..2000 {
        let (level, l0) = (rng.uniform(40.0, 140.0), rng.uniform(80.0, 90.0));
        let (q, tc, t) = (
            rng.uniform(2.0, 6.0),
            rng.uniform(1.0, 40.0),
            rng.uniform(0.1, 10.0),
        );
        let half = allowed_time(level, l0, q, tc).unwrap() / 2.0;
        let next = allowed_time(level + q, l0, q, tc).unwrap();
        assert!((next - half).abs() <= 1e-13 * half, "{level} {l0} {q}");
        let dose = noise_dose(&[level], &[t], l0, q, tc).unwrap();
        let same = noise_dose(&[level + q], &[t / 2.0], l0, q, tc).unwrap();
        assert!((same - dose).abs() <= 1e-13 * dose);
        let full_day = noise_dose(&[level], &[tc], l0, q, tc).unwrap();
        assert!((twa_from_dose(full_day, l0, q).unwrap() - level).abs() < 1e-12);
    }
}

/// With an exchange rate of exactly 10·log10 2 dB (equal energy), the dose of
/// mixed segments equals the dose of their Leq over the same time. The 3 dB
/// rule is close: for this 8-hour day (3 h at 82, 0.5 h at 97, 1.5 h at 90 and
/// 3 h at 75 dBA) its TWA is 87.592 dBA against an Leq of 87.574 dBA.
#[test]
fn equal_energy_exchange_rate_matches_leq() {
    let q = power_to_db(2.0).unwrap();
    let (lv, hours) = ([82.0, 97.0, 90.0, 75.0], [3.0, 0.5, 1.5, 3.0]);
    let level = leq(&lv, &hours).unwrap();
    let exact = noise_dose(&lv, &hours, 85.0, q, 8.0).unwrap();
    let steady = noise_dose(&[level], &[8.0], 85.0, q, 8.0).unwrap();
    assert!((exact - steady).abs() < 1e-12 * exact);
    let twa = twa_from_dose(noise_dose(&lv, &hours, 85.0, 3.0, 8.0).unwrap(), 85.0, 3.0).unwrap();
    assert_eq!(
        [(twa * 1000.0).round(), (level * 1000.0).round()],
        [87_592.0, 87_574.0]
    );
}

#[test]
fn normalize_reaches_the_target_or_the_ceiling() {
    let sr = 8000;
    let tone = osc(220.0, 0.25, sr, 0.0).unwrap();
    for (gain, target, ceiling) in [(0.01, -20.0, -1.0), (0.9, -3.0, -1.0), (0.3, -12.0, -6.0)] {
        let x: Vec<f64> = tone.iter().map(|v| gain * v).collect();
        let y = normalize(&x, 1, target, ceiling, sr).unwrap();
        let (rms, peak) = (rms_dbfs(&y).unwrap(), peak_dbfs(&y).unwrap());
        assert!(peak <= ceiling + 1e-12);
        assert!((rms - target).abs() < 1e-9 || (peak - ceiling).abs() < 1e-9);
    }
}

/// The loudness meter's scaling laws: a gain of g adds 20·log10 g LU, a copy of
/// a mono channel in a second channel adds 10·log10 2, and the result does not
/// depend on the order of the channels.
#[test]
fn loudness_scales_with_gain_and_channels() {
    let sr = 44_100;
    let x: Vec<f64> = osc(440.0, 1.5, sr, 0.0)
        .unwrap()
        .iter()
        .zip(osc(97.0, 1.5, sr, 0.3).unwrap())
        .map(|(a, b)| 0.2 * a + 0.1 * b)
        .collect();
    let mono = integrated_loudness(&x, &[1.0], sr).unwrap();
    for g in [0.01, 0.5, 3.0] {
        let y: Vec<f64> = x.iter().map(|v| g * v).collect();
        let louder = integrated_loudness(&y, &[1.0], sr).unwrap();
        assert!((louder - mono - amplitude_to_db(g).unwrap()).abs() < 1e-10);
    }
    let twice: Vec<f64> = x.iter().flat_map(|&v| [v, v]).collect();
    let stereo = integrated_loudness(&twice, &[1.0, 1.0], sr).unwrap();
    assert!((stereo - mono - power_to_db(2.0).unwrap()).abs() < 1e-10);
    let quiet: Vec<f64> = x.iter().flat_map(|&v| [v, 0.25 * v]).collect();
    let swapped: Vec<f64> = x.iter().flat_map(|&v| [0.25 * v, v]).collect();
    let (a, b) = (
        integrated_loudness(&quiet, &[1.0, SURROUND_WEIGHT], sr).unwrap(),
        integrated_loudness(&swapped, &[SURROUND_WEIGHT, 1.0], sr).unwrap(),
    );
    assert!((a - b).abs() < 1e-12);
}

/// |H(e^{iω})| of a biquad `[b0, b1, b2, a1, a2]` at `hz`, in dB.
fn response_db([b0, b1, b2, a1, a2]: [f64; 5], hz: f64, sr: u32) -> f64 {
    let w = 2.0 * PI * hz / f64::from(sr);
    let (c1, s1, c2, s2) = (w.cos(), w.sin(), (2.0 * w).cos(), (2.0 * w).sin());
    let num = (b0 + b1 * c1 + b2 * c2).hypot(b1 * s1 + b2 * s2);
    let den = (1.0 + a1 * c1 + a2 * c2).hypot(a1 * s1 + a2 * s2);
    20.0 * (num / den).log10()
}

/// The K-weighting is the same analogue filter at every rate: below 6 kHz its
/// response at 44.1, 96 and 192 kHz stays within 0.01 dB of the 48 kHz filter
/// (the residue is the bilinear transform's frequency warping). At 1 kHz it
/// adds about +0.70 dB, so a 1 kHz sine at -23 dBFS in two channels reads
/// -22.99 LUFS.
#[test]
fn k_weighting_has_the_same_response_at_every_rate() {
    let total = |sr: u32, hz: f64| -> f64 {
        let [shelf, hp] = k_weighting(sr).unwrap();
        response_db(shelf, hz, sr) + response_db(hp, hz, sr)
    };
    let mut worst: f64 = 0.0;
    for sr in [44_100, 96_000, 192_000] {
        for hz in (1..=600).map(|k| 10.0 * f64::from(k)) {
            worst = worst.max((total(sr, hz) - total(48_000, hz)).abs());
        }
    }
    println!("largest K-weighting difference from 48 kHz below 6 kHz: {worst:.5} dB");
    assert!(worst < 0.01);
    let at_1k = total(48_000, 1000.0);
    assert!((at_1k - 0.6977).abs() < 0.001, "{at_1k}");
}

#[test]
fn nan_and_infinity_are_rejected() {
    let bad = [f64::NAN, f64::INFINITY, f64::NEG_INFINITY];
    for v in bad {
        assert!(amplitude_to_db(v).is_err() && power_to_db(v).is_err());
        assert!(db_to_amplitude(v).is_err() && db_to_power(v).is_err());
        assert!(sum_incoherent_db(&[60.0, v]).is_err());
        assert!(energetic_mean_db(&[v]).is_err());
        assert!(leq(&[60.0], &[v]).is_err() && leq(&[v], &[1.0]).is_err());
        assert!(rms_dbfs(&[0.1, v]).is_err() && peak_dbfs(&[v]).is_err());
        assert!(normalize(&[0.1, v], 1, -20.0, -1.0, 8000).is_err());
        assert!(normalize(&[0.1, 0.2], 1, v, -1.0, 8000).is_err());
        assert!(allowed_time(v, 85.0, 3.0, 8.0).is_err());
        assert!(allowed_time(85.0, v, 3.0, 8.0).is_err());
        assert!(noise_dose(&[v], &[1.0], 85.0, 3.0, 8.0).is_err());
        assert!(twa_from_dose(v, 85.0, 3.0).is_err());
        assert!(k_weight(&[0.1, v], 1, 48_000).is_err());
        assert!(integrated_loudness(&[v; 2000], &[1.0], 4000).is_err());
        assert!(integrated_loudness(&[0.1; 2000], &[v], 4000).is_err());
    }
    // Finite input whose energy overflows an f64 is an error too, not NaN or ∞.
    assert!(integrated_loudness(&[1e200; 2000], &[1.0], 4000).is_err());
    assert!(short_term_loudness(&[1e200; 12_000], &[1.0], 4000).is_err());
    assert!(normalize(&[0.3; 1000], 1, -20.0, -1.0, 1000).is_err()); // a constant leaves only rounding
    // A silent window is -inf in the meters, a number in no vector file.
    let gap: Vec<f64> = [vec![0.0; 2000], vec![0.3; 2000]].concat();
    let m = momentary_loudness(&gap, &[1.0], 4000).unwrap();
    assert_eq!(m[0], f64::NEG_INFINITY);
    assert!(m.last().unwrap().is_finite());
}
