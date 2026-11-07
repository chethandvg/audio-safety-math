// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Loudness per ITU-R BS.1770: K-weighting, gating blocks and integrated,
//! momentary and short-term loudness in LUFS.

use crate::{Result, SafetyError, all_finite, each, frames, require};

/// The absolute gate of BS.1770: blocks at or below −70 LKFS are ignored.
pub const ABSOLUTE_GATE_LKFS: f64 = -70.0;

/// The relative gate of BS.1770: blocks at or below 10 LU under the loudness of
/// the blocks that pass the absolute gate are ignored.
pub const RELATIVE_GATE_LU: f64 = -10.0;

/// The BS.1770 channel weight of the left and right surround channels (+1.5 dB);
/// left, right and centre have weight 1, and the LFE channel is left out (weight 0).
pub const SURROUND_WEIGHT: f64 = 1.41;

/// The lowest sample rate the K-weighting accepts: above twice the shelf's
/// 1681.97 Hz, so that the shelf lies below the Nyquist frequency.
pub const K_WEIGHTING_MIN_SR: u32 = 3364;

/// `−0.691 + 10 · log₁₀ p`: the loudness of a weighted mean square `p`
/// (−∞ for 0). The offset makes a 997 Hz sine read the same in LKFS as its RMS
/// level in each of two channels.
const OFFSET_LU: f64 = -0.691;

/// Hops of 100 ms: a 400 ms block (75 % overlap) is 4 hops, a 3 s window 30.
const HOPS_PER_SECOND: u64 = 10;
const BLOCK_HOPS: usize = 4;
const SHORT_TERM_HOPS: usize = 30;

/// A second-order section `H(z) = (b₀ + b₁z⁻¹ + b₂z⁻²) / (1 + a₁z⁻¹ + a₂z⁻²)`,
/// stored as `[b₀, b₁, b₂, a₁, a₂]`.
type Biquad = [f64; 5];

/// The analogue prototype of one K-weighting stage, in the frequency `s / ω₀`
/// normalised to its corner `f0`:
/// `H(s) = (vh·s² + (vb/q)·s + vl) / (s² + s/q + 1)`.
struct Prototype {
    f0: f64,
    q: f64,
    vl: f64,
    vb: f64,
    vh: f64,
}

/// The two stages of BS.1770's K-weighting as analogue prototypes. They are
/// the published 48 kHz coefficients mapped back through the bilinear
/// transform (vectors/generate.py derives them in 50-digit arithmetic), so
/// [`k_weighting`]`(48000)` returns those coefficients again: a high shelf of
/// +4.0 dB above about 1.7 kHz (the acoustic effect of the head), then a
/// high-pass at about 38 Hz (the revised low-frequency B-weighting curve).
const SHELF: Prototype = Prototype {
    f0: 1681.9744509555328,
    q: 0.7071752369554196,
    vl: 1.0,
    vb: 1.2587209302325617,
    vh: 1.5848647011308554,
};
const HIGH_PASS: Prototype = Prototype {
    f0: 38.13547087602444,
    q: 0.5003270373238773,
    vl: 0.0,
    vb: 0.0,
    vh: 1.0049948987146884,
};

impl Prototype {
    /// The bilinear transform at `sr`, prewarped so that `f0` maps exactly:
    /// `s ↦ (1/K)·(1 − z⁻¹)/(1 + z⁻¹)` with `K = tan(π·f0/sr)`.
    fn biquad(&self, sr: f64) -> Biquad {
        let k = (std::f64::consts::PI * self.f0 / sr).tan();
        let (kq, k2) = (k / self.q, k * k);
        let d = 1.0 + kq + k2;
        [
            (self.vh + self.vb * kq + self.vl * k2) / d,
            2.0 * (self.vl * k2 - self.vh) / d,
            (self.vh - self.vb * kq + self.vl * k2) / d,
            2.0 * (k2 - 1.0) / d,
            (1.0 - kq + k2) / d,
        ]
    }
}

fn loudness_rate(sr: u32) -> Result<f64> {
    require(
        sr >= K_WEIGHTING_MIN_SR,
        "sr",
        f64::from(sr),
        "at least K_WEIGHTING_MIN_SR (3364) samples per second",
    )
    .map(|()| f64::from(sr))
}

/// The K-weighting filter of ITU-R BS.1770 at sample rate `sr`: the high shelf
/// and the high-pass, each as `[b₀, b₁, b₂, a₁, a₂]` with `a₀ = 1`. At 48 kHz
/// these are the standard's published coefficients (to within 10⁻¹⁵); at other
/// rates they are the same analogue filters, discretised by the prewarped
/// bilinear transform. `sr` must be at least [`K_WEIGHTING_MIN_SR`].
pub fn k_weighting(sr: u32) -> Result<[[f64; 5]; 2]> {
    let fs = loudness_rate(sr)?;
    Ok([SHELF.biquad(fs), HIGH_PASS.biquad(fs)])
}

/// Filters channel `c` of the interleaved `y` in place (transposed direct form II).
fn filter_channel(y: &mut [f64], channels: usize, c: usize, [b0, b1, b2, a1, a2]: Biquad) {
    let (mut s1, mut s2) = (0.0, 0.0);
    for v in y.iter_mut().skip(c).step_by(channels) {
        let x = *v;
        let out = b0 * x + s1;
        s1 = b1 * x - a1 * out + s2;
        s2 = b2 * x - a2 * out;
        *v = out;
    }
}

/// `x` through both stages, each channel from rest; the output is not checked.
fn k_filtered(x: &[f64], channels: usize, sr: u32) -> Result<Vec<f64>> {
    frames("x", x, channels)?;
    all_finite("x", x)?;
    let stages = k_weighting(sr)?;
    let mut y = x.to_vec();
    for c in 0..channels {
        stages
            .iter()
            .for_each(|&s| filter_channel(&mut y, channels, c, s));
    }
    Ok(y)
}

/// `x` (interleaved, `channels` channels) through the [`k_weighting`] filter,
/// each channel separately, starting from rest. An output sample too large for
/// an `f64` (input near the largest `f64`) is an error.
pub fn k_weight(x: &[f64], channels: usize, sr: u32) -> Result<Vec<f64>> {
    let y = k_filtered(x, channels, sr)?;
    let fit = "small enough to fit in an f64 in every element";
    each("result", &y, f64::is_finite, fit).map(|()| y)
}

/// The requirement on a signal whose loudness is asked for.
const ENERGY: &str = "small enough that its K-weighted energy, weighted by channel, fits in an f64";

/// A mean square that must fit in an `f64`.
fn energy(p: f64) -> Result<f64> {
    require(p.is_finite(), "x", p, ENERGY).map(|()| p)
}

/// The K-weighted energy of `x` in 100 ms hops: `energy[k·channels + c]` is the
/// sum of squares of channel `c` over frames `bounds[k] .. bounds[k + 1]`, with
/// `bounds[k] = ⌊k·sr/10⌋`. Only complete hops are kept.
struct Hops {
    energy: Vec<f64>,
    bounds: Vec<usize>,
    channels: usize,
}

impl Hops {
    fn new(x: &[f64], weights: &[f64], sr: u32) -> Result<Self> {
        let channels = weights.len();
        require(channels >= 1, "weights", 0.0, "at least one channel weight")?;
        each(
            "weights",
            weights,
            |w| w.is_finite() && w >= 0.0,
            "finite and at least 0 for every channel",
        )?;
        let n = frames("x", x, channels)? as u64;
        let y = k_filtered(x, channels, sr)?;
        let hops = n * HOPS_PER_SECOND / u64::from(sr);
        let bounds: Vec<usize> = (0..=hops)
            .map(|k| (k * u64::from(sr) / HOPS_PER_SECOND) as usize)
            .collect();
        let mut energy = vec![0.0; hops as usize * channels];
        for (k, hop) in bounds.windows(2).enumerate() {
            let sums = &mut energy[k * channels..(k + 1) * channels];
            for frame in y[hop[0] * channels..hop[1] * channels].chunks_exact(channels) {
                sums.iter_mut().zip(frame).for_each(|(s, v)| *s += v * v);
            }
        }
        Ok(Self {
            energy,
            bounds,
            channels,
        })
    }

    /// The channel-weighted mean square `Σ_c G_c · z_c` of every window of
    /// `hops` consecutive hops, one window per hop. Channels of weight 0 take
    /// no part; a mean square too large for an `f64` is an error.
    fn window_powers(&self, weights: &[f64], hops: usize) -> Result<Vec<f64>> {
        let count = (self.bounds.len() - 1).saturating_sub(hops - 1);
        let powers: Vec<f64> = (0..count)
            .map(|j| {
                let len = (self.bounds[j + hops] - self.bounds[j]) as f64;
                let rows = &self.energy[j * self.channels..(j + hops) * self.channels];
                let weighted: f64 = weights
                    .iter()
                    .enumerate()
                    .filter(|&(_, &g)| g > 0.0)
                    .map(|(c, g)| g * rows.iter().skip(c).step_by(self.channels).sum::<f64>())
                    .sum();
                weighted / len
            })
            .collect();
        each("x", &powers, f64::is_finite, ENERGY).map(|()| powers)
    }
}

fn lufs(power: f64) -> f64 {
    OFFSET_LU + 10.0 * power.log10()
}

fn mean(xs: impl ExactSizeIterator<Item = f64>) -> f64 {
    let n = xs.len() as f64;
    xs.sum::<f64>() / n
}

/// The gated integrated loudness of ITU-R BS.1770 in LUFS. `x` holds
/// `weights.len()` interleaved channels and `weights[c]` is the weight of
/// channel `c` (1 for left, right and centre, [`SURROUND_WEIGHT`] for the
/// surround channels, 0 for LFE). Each channel is K-weighted; the mean squares
/// `z` over blocks of 400 ms, one every 100 ms (75 % overlap), give block
/// loudnesses `l = −0.691 + 10·log₁₀ Σ G_c·z_c`; the blocks above
/// [`ABSOLUTE_GATE_LKFS`] give a relative gate [`RELATIVE_GATE_LU`] below
/// their own loudness, and the result is the loudness of the mean power of the
/// blocks above both gates. A signal shorter than 400 ms, one whose blocks all
/// fall below the absolute gate (silence, say), or one whose K-weighted energy
/// overflows an `f64` (samples of about 10¹⁵² or more) is an error.
pub fn integrated_loudness(x: &[f64], weights: &[f64], sr: u32) -> Result<f64> {
    let blocks = Hops::new(x, weights, sr)?.window_powers(weights, BLOCK_HOPS)?;
    require(
        !blocks.is_empty(),
        "x",
        x.len() as f64,
        "at least 400 ms long (one gating block)",
    )?;
    let level: Vec<f64> = blocks.iter().map(|&p| lufs(p)).collect();
    let above = |gate: f64| {
        let kept = blocks.iter().zip(&level).filter(move |&(_, &l)| l > gate);
        kept.map(|(&p, _)| p).collect::<Vec<f64>>()
    };
    let absolute = above(ABSOLUTE_GATE_LKFS);
    if absolute.is_empty() {
        return Err(SafetyError {
            name: "x",
            value: level.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            requirement: "louder than the absolute gate (-70 LKFS) in at least one 400 ms block",
        });
    }
    let relative = lufs(energy(mean(absolute.into_iter()))?) + RELATIVE_GATE_LU;
    let gated = above(relative.max(ABSOLUTE_GATE_LKFS));
    Ok(lufs(energy(mean(gated.into_iter()))?))
}

/// Momentary loudness in LUFS: the loudness of every 400 ms window, one every
/// 100 ms (the first covers 0–0.4 s), ungated; a silent window is −∞, and so is
/// one whose squares all underflow to 0 (samples below about 10⁻¹⁶²). Channels and
/// weights are as in [`integrated_loudness`], and so is the error for energy
/// that overflows; a signal shorter than 400 ms gives no values.
pub fn momentary_loudness(x: &[f64], weights: &[f64], sr: u32) -> Result<Vec<f64>> {
    meter(x, weights, sr, BLOCK_HOPS)
}

/// Short-term loudness in LUFS: the loudness of every 3 s window, one every
/// 100 ms, ungated; otherwise as [`momentary_loudness`].
pub fn short_term_loudness(x: &[f64], weights: &[f64], sr: u32) -> Result<Vec<f64>> {
    meter(x, weights, sr, SHORT_TERM_HOPS)
}

fn meter(x: &[f64], weights: &[f64], sr: u32, hops: usize) -> Result<Vec<f64>> {
    let windows = Hops::new(x, weights, sr)?.window_powers(weights, hops)?;
    Ok(windows.into_iter().map(lufs).collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// ITU-R BS.1770, tables 1 and 2: the 48 kHz coefficients `[b₀, b₁, b₂, a₁, a₂]`.
    const PUBLISHED: [[f64; 5]; 2] = [
        [
            1.53512485958697,
            -2.69169618940638,
            1.19839281085285,
            -1.69065929318241,
            0.73248077421585,
        ],
        [1.0, -2.0, 1.0, -1.99004745483398, 0.99007225036621],
    ];

    #[test]
    fn k_weighting_at_48_khz_is_the_published_filter() {
        let got = k_weighting(48_000).unwrap();
        for (g, p) in got.iter().flatten().zip(PUBLISHED.iter().flatten()) {
            assert!((g - p).abs() < 1e-14, "{g} vs {p}");
        }
    }

    #[test]
    fn k_weighting_keeps_dc_out_and_the_shelf_gain() {
        for sr in [K_WEIGHTING_MIN_SR, 8000, 44_100, 48_000, 96_000, 192_000] {
            let [shelf, hp] = k_weighting(sr).unwrap();
            let at = |[b0, b1, b2, a1, a2]: Biquad, z: f64| {
                (b0 + b1 * z + b2 * z * z) / (1.0 + a1 * z + a2 * z * z)
            };
            assert!(
                (at(shelf, 1.0) - 1.0).abs() < 1e-12,
                "shelf DC gain at {sr}"
            );
            assert!(at(hp, 1.0).abs() < 1e-15, "high-pass DC gain at {sr}");
            // Just above the minimum rate the shelf's corner is at the Nyquist
            // frequency, and the gain there loses digits to cancellation.
            let tol = if sr < 8000 { 1e-6 } else { 1e-12 };
            assert!(
                (at(shelf, -1.0) - SHELF.vh).abs() < tol,
                "shelf Nyquist gain at {sr}: {}",
                at(shelf, -1.0) - SHELF.vh
            );
        }
        assert!(k_weighting(K_WEIGHTING_MIN_SR - 1).is_err());
    }

    #[test]
    fn short_and_silent_signals() {
        let sr = 8000;
        assert!(integrated_loudness(&[0.5; 3199], &[1.0], sr).is_err()); // under 400 ms
        assert!(integrated_loudness(&[0.0; 8000], &[1.0], sr).is_err());
        let m = momentary_loudness(&[0.0; 4000], &[1.0], sr).unwrap();
        assert_eq!(m, [f64::NEG_INFINITY, f64::NEG_INFINITY]);
        assert!(
            short_term_loudness(&[0.5; 8000], &[1.0], sr)
                .unwrap()
                .is_empty()
        );
        assert!(integrated_loudness(&[0.5; 8000], &[], sr).is_err());
        assert!(integrated_loudness(&[0.5; 8000], &[-1.0], sr).is_err());
        assert!(integrated_loudness(&[0.5; 8001], &[1.0, 1.0], sr).is_err());
        assert!(k_weight(&[f64::NAN], 1, sr).is_err());
    }

    /// Finite input whose filtered output or energy overflows is an error, never NaN or ∞.
    #[test]
    fn overflowing_energy_is_an_error() {
        let huge = [1e200; 1600];
        assert_eq!(
            integrated_loudness(&huge, &[1.0], 4000).unwrap_err().name(),
            "x"
        );
        assert!(momentary_loudness(&huge, &[1.0], 4000).is_err());
        let tone: Vec<f64> = (0..3200)
            .map(|i| 0.3 * (0.7 * f64::from(i)).sin())
            .collect();
        assert!(integrated_loudness(&tone, &[1e308, 1e308], 4000).is_err());
        let edge: Vec<f64> = [1e308, -1e308].repeat(4);
        assert_eq!(k_weight(&edge, 1, 48_000).unwrap_err().name(), "result");
        // A channel of weight 0 takes no part, however loud it is.
        let lfe: Vec<f64> = tone.iter().flat_map(|&v| [v, 1e200]).collect();
        let mono = integrated_loudness(&tone, &[1.0], 4000).unwrap();
        assert_eq!(integrated_loudness(&lfe, &[1.0, 0.0], 4000).unwrap(), mono);
    }
}
