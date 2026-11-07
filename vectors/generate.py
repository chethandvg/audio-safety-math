# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
"""Write the golden vectors for audio-safety-math from independent references.

Nothing here calls the Rust code. The references are:

* decibels, sums of levels, Leq and the dose rule: the textbook formulas,
  evaluated directly in 50-digit arithmetic with mpmath (without the crate's
  factoring around the largest term);
* rms_dbfs and normalize: the functions of the project's earlier toolkit,
  lifesound.py, copied below with its box_lowpass and highpass, with two
  deliberate changes: rms_dbfs no longer adds 1e-12 before the logarithm, and
  normalize passes its sample rate to highpass (rms_dbfs of samples whose
  squares underflow or overflow, and peak_dbfs, in mpmath);
* the K-weighting coefficients: the 48 kHz coefficients published in ITU-R
  BS.1770 (tables 1 and 2), mapped back through the bilinear transform to two
  analogue prototypes and discretised again at each sample rate, in 50-digit
  arithmetic;
* k_weight: scipy.signal.lfilter with those coefficients;
* integrated, momentary and short-term loudness: pyloudnorm's meter (its
  blocks and gates) with those coefficients in place of its own filters,
  cross-checked against `bs1770`, a direct transcription of the BS.1770
  formulas, which is also the reference for one case pyloudnorm cannot express
  (an LFE channel of weight 0).

It also writes tests/ebu3341.json: the EBU Tech 3341 test signals used by
tests/ebu3341.rs, with the published expected values and pyloudnorm's
readings of the same signals (its "DeMan" filters, which match BS.1770 at
48 kHz). Run with numpy, scipy, mpmath and pyloudnorm installed:

    python vectors/generate.py
"""

from __future__ import annotations

import json
from pathlib import Path

import mpmath as mp
import numpy as np
import pyloudnorm
import scipy.signal

mp.mp.dps = 50
HERE = Path(__file__).resolve().parent
TESTS = HERE.parent / "tests"
SR = 44100
CLOSE = {"rel": 1e-12, "abs": 1e-12}
RATIO = {"rel": 1e-12, "abs": 0}
LU = {"rel": 0, "abs": 1e-9}


def write(name: str, description: str, cases: list[dict], tolerance: dict = CLOSE) -> None:
    """One case per line, so that diffs show which cases changed."""
    doc = {"schema": "golden-vectors@1", "description": description, "tolerance": tolerance, "cases": cases}
    head = json.dumps({k: doc[k] for k in ("schema", "description", "tolerance")}, indent=1, ensure_ascii=False)
    body = ",\n".join("  " + json.dumps(c, ensure_ascii=False) for c in cases)
    text = head[:-2] + ',\n "cases": [\n' + body + "\n ]\n}\n"
    assert json.loads(text) == doc
    (HERE / f"{name}.json").write_text(text, encoding="utf-8", newline="\n")  # LF on every platform


def case(claim: str, fn: str, args: list, expected, note: str | None = None, tolerance: dict | None = None) -> dict:
    c = {"claim": claim, "fn": fn, "args": args, "expected": expected}
    if tolerance is not None:
        c["tolerance"] = tolerance
    if note:
        c["note"] = note
    return c


def error(claim: str, fn: str, args: list, note: str) -> dict:
    return {"claim": claim, "fn": fn, "args": args, "error": True, "note": note}


def f(x) -> float:
    return float(x)


# Decibels and levels, in 50-digit arithmetic ----------------------------------------


def db(ratio, scale) -> float:
    return f(scale * mp.log10(mp.mpf(ratio)))


def undb(level, scale) -> float:
    return f(mp.power(10, mp.mpf(level) / scale))


def power_sum_db(levels, weights=None) -> float:
    weights = weights or [1] * len(levels)
    total = mp.fsum(mp.mpf(w) * mp.power(10, mp.mpf(lv) / 10) for lv, w in zip(levels, weights, strict=True))
    return f(10 * mp.log10(total))


def leq(levels, durations) -> float:
    total = mp.fsum(mp.mpf(t) for t in durations)
    return power_sum_db(levels, [mp.mpf(t) / total for t in durations])


def decibel_vectors() -> None:
    db_claim, sum_claim = "safety:DB-CONVERSIONS", "safety:LEVEL-SUMS"
    cases = [
        case(db_claim, "amplitude_to_db", [r], db(r, 20)) for r in (1.0, 2.0, 10.0, 0.5, 0.001, 3.7, 1e-300, 1e300)
    ]
    cases += [case(db_claim, "power_to_db", [r], db(r, 10)) for r in (1.0, 2.0, 100.0, 0.25, 1e-12, 1.5e200)]
    cases += [
        case(db_claim, "db_to_amplitude", [x], undb(x, 20), tolerance=RATIO)
        for x in (0.0, 6.0, -6.0, 20.0, -120.0, 3.5)
    ]
    cases += [
        case(db_claim, "db_to_power", [x], undb(x, 10), tolerance=RATIO) for x in (0.0, 3.0, -3.0, 10.0, 90.0, -250.5)
    ]
    sums = ([80.0, 80.0], [90.0, 80.0, 70.0], [60.0] * 10, [-20.0, 400.0], [400.0, 400.0], [55.5], [-3000.0, -3001.0])
    cases += [case(sum_claim, "sum_incoherent_db", [lv], power_sum_db(lv)) for lv in sums]
    for lv in ([80.0, 80.0], [90.0, 80.0, 70.0], [60.0, 70.0], [1000.0, 990.0, 1005.0], [-12.5]):
        n = len(lv)
        cases.append(case(sum_claim, "energetic_mean_db", [lv], power_sum_db(lv, [mp.mpf(1) / n] * n)))
    for lv, t in (
        ([85.0, 91.0, 70.0], [4.0, 2.0, 2.0]),
        ([94.0], [0.5]),
        ([100.0, 60.0], [0.25, 7.75]),
        ([75.0, 120.0, 80.0], [3.0, 0.0, 5.0]),
        ([60.0, 63.0], [1.0, 1.0]),
        ([60.0, 70.0], [1e308, 1e308]),  # the durations' sum exceeds the f64 range
        ([0.0, 6000.0], [1e10, 1e-320]),  # 1e-330 of the longest duration, at a level that outweighs it
    ):
        cases.append(case(sum_claim, "leq", [lv, t], leq(lv, t)))
    write(
        "decibels", "Ratios to and from decibels, incoherent sums, energetic means and Leq (mpmath, 50 digits).", cases
    )


# Noise dose, in 50-digit arithmetic ----------------------------------------------


def allowed(level, l0, q, tc):
    return mp.mpf(tc) / mp.power(2, (mp.mpf(level) - l0) / q)


def dose_vectors() -> None:
    claim = "safety:DOSE-MATCHES-REFERENCE"
    niosh, osha = (85.0, 3.0, 8.0), (90.0, 5.0, 8.0)
    cases = [
        case(claim, "allowed_time", [lv, *rule], f(allowed(lv, *rule)), note=note)
        for lv, rule, note in (
            (85.0, niosh, "NIOSH: 8 h at the criterion level"),
            (88.0, niosh, "4 h"),
            (94.0, niosh, "1 h"),
            (100.0, niosh, "15 min"),
            (82.0, niosh, "16 h"),
            (87.3, niosh, None),
            (90.0, osha, "OSHA parameters: 8 h"),
            (95.0, osha, "4 h"),
            (105.0, osha, "1 h"),
            (92.5, osha, None),
            (80.0, (80.0, 3.0, 40.0), "a 40-hour weekly allowance at 80 dB"),
            (90.0, (80.0, 3.0, 40.0), None),
        )
    ]
    for lv, hours, rule in (
        ([88.0, 85.0], [2.0, 4.0], niosh),
        ([95.0, 90.0, 85.0, 80.0], [1.0, 2.0, 3.0, 2.0], niosh),
        ([100.0], [0.25], niosh),
        ([92.0, 97.5], [3.0, 0.5], osha),
        ([84.2, 101.3, 70.0], [5.5, 0.1, 2.4], niosh),
        ([], [], niosh),
        ([5000.0, 85.0], [0.0, 8.0], niosh),  # no time at 5000 dBA is no dose
        ([3160.0], [1e-10], niosh),  # 2^(3075/3) overflows an f64; the dose does not
    ):
        total = mp.fsum(mp.mpf(c) / allowed(lv_i, *rule) for lv_i, c in zip(lv, hours, strict=True))
        cases.append(case(claim, "noise_dose", [lv, hours, *rule], f(100 * total)))
    for dose, rule in (
        (100.0, niosh),
        (200.0, niosh),
        (50.0, niosh),
        (400.0, niosh),
        (137.0, niosh),
        (200.0, osha),
        (3.2, osha),
    ):
        l0, q, _ = rule
        cases.append(case(claim, "twa_from_dose", [dose, l0, q], f(l0 + q * mp.log(mp.mpf(dose) / 100, 2))))
    write("dose", "Allowed time, daily dose and TWA by the exchange-rate rule (mpmath, 50 digits).", cases)


# The reference implementation: lifesound.py ------------------------------------------


def box_lowpass(x: np.ndarray, cutoff_hz: float, sr: int = SR, passes: int = 3) -> np.ndarray:
    h = max(int(0.44 * sr / cutoff_hz) // 2, 1)
    L = 2 * h + 1  # noqa: N806 - copied unchanged
    y = x.astype(float)
    for _ in range(passes):
        padded = np.pad(y, [(h + 1, h)] + [(0, 0)] * (y.ndim - 1), mode="edge")
        c = np.cumsum(padded, axis=0)
        y = (c[L:] - c[:-L]) / L
    return y


def highpass(x: np.ndarray, cutoff_hz: float = 12.0, sr: int = SR) -> np.ndarray:
    return x - box_lowpass(x, cutoff_hz, sr)


def rms_dbfs(x: np.ndarray) -> float:
    return 20 * np.log10(np.sqrt(np.mean(np.square(x))))  # lifesound.py: ... + 1e-12


def normalize(x: np.ndarray, target_rms_dbfs: float = -20.0, peak_ceiling_dbfs: float = -1.0, sr: int = SR):
    x = highpass(x, 12.0, sr)  # lifesound.py: highpass(x, 12.0), at its default sr
    y = x * 10 ** ((target_rms_dbfs - rms_dbfs(x)) / 20)
    peak = np.max(np.abs(y))
    ceiling = 10 ** (peak_ceiling_dbfs / 20)
    if peak > ceiling:
        y *= ceiling / peak
    return y


def tones(sr: int, seconds: float, parts, decimals: int = 4) -> np.ndarray:
    """A sum of sines amp * sin(2 pi f t + phase), rounded to `decimals` places to keep the files small."""
    t = np.arange(round(seconds * sr)) / sr
    return np.round(sum(a * np.sin(2 * np.pi * hz * t + ph) for a, hz, ph in parts), decimals)


def interleave(*channels: np.ndarray) -> np.ndarray:
    return np.stack(channels, axis=1).reshape(-1)


def level_vectors() -> None:
    claim = "safety:LEVELS-MATCH-REFERENCE"
    mono = np.round(tones(1000, 0.2, [(0.5, 50.0, 0.3), (0.2, 130.0, 1.0)]) + 0.1, 4)  # DC for the high-pass to remove
    quiet = tones(1000, 0.15, [(0.01, 70.0, 0.0), (0.004, 210.0, 2.0)], 6)
    left, right = tones(1000, 0.1, [(0.6, 40.0, 0.0)]), tones(1000, 0.1, [(0.3, 90.0, 0.5), (0.2, 300.0, 1.5)])
    stereo = interleave(left, right)
    spiky = np.zeros(120)
    spiky[[10, 50, 51]] = [0.9, -0.4, 0.3]
    cases = []
    for name, x in (("mono", mono), ("quiet", quiet), ("stereo", stereo), ("spiky", spiky)):
        cases.append(case(claim, "rms_dbfs", [x.tolist()], f(rms_dbfs(x)), note=name))
        peak = mp.log10(max(abs(mp.mpf(v)) for v in x.tolist()))
        cases.append(case(claim, "peak_dbfs", [x.tolist()], f(20 * peak), note=name))
    # Samples whose squares underflow or overflow an f64 (numpy's rms_dbfs fails here): reference mpmath.
    for x, note in (([1e-200, -3e-200, 2e-200], "squares underflow"), ([1e300, -1e300, 5e299], "squares overflow")):
        mean_square = mp.fsum(mp.mpf(v) ** 2 for v in x) / len(x)
        cases.append(case(claim, "rms_dbfs", [x], f(10 * mp.log10(mean_square)), note=note + "; reference: mpmath"))
    for x, channels, target, ceiling, sr, note in (
        (mono, 1, -20.0, -1.0, 1000, "the RMS target is reached"),
        (mono, 1, -6.0, -3.0, 1000, "the peak ceiling binds"),
        (quiet, 1, -20.0, -1.0, 1000, "a quiet signal is raised"),
        (stereo, 2, -18.0, -1.0, 1000, "stereo: one RMS over both channels"),
        (spiky, 1, -20.0, -1.0, 1000, "a few clicks: the ceiling binds"),
        (mono, 1, -20.0, -1.0, 8000, "the same samples at 8 kHz: a longer high-pass"),
    ):
        frames = x.reshape(-1, channels) if channels > 1 else x
        y = normalize(frames, target, ceiling, sr).reshape(-1)
        cases.append(case(claim, "normalize", [x.tolist(), channels, target, ceiling, sr], y.tolist(), note=note))
    write("levels", "rms_dbfs, peak_dbfs and normalize: lifesound.py's functions without the 1e-12 (numpy).", cases)


# K-weighting: BS.1770's published filter and its analogue prototypes ---------------------

# ITU-R BS.1770, tables 1 and 2: the K-weighting stages at 48 kHz as [b0, b1, b2, a1, a2], a0 = 1.
PUBLISHED = {
    "shelf": ["1.53512485958697", "-2.69169618940638", "1.19839281085285", "-1.69065929318241", "0.73248077421585"],
    "high_pass": ["1.0", "-2.0", "1.0", "-1.99004745483398", "0.99007225036621"],
}


def prototype(coefficients: list[str], fs: int = 48000) -> dict:
    """Invert the prewarped bilinear transform of H(s) = (vh s^2 + (vb/q) s + vl) / (s^2 + s/q + 1)."""
    b0, b1, b2, a1, a2 = (mp.mpf(c) for c in coefficients)
    k2 = (1 + a1 + a2) / (1 - a1 + a2)
    d = 4 / (1 - a1 + a2)
    kq = d * (1 - a2) / 2
    k = mp.sqrt(k2)
    return {
        "f0": fs / mp.pi * mp.atan(k),
        "q": k / kq,
        "vl": d * (b0 + b1 + b2) / (4 * k2),
        "vb": d * (b0 - b2) / (2 * kq),
        "vh": d * (b0 - b1 + b2) / 4,
    }


PROTOTYPES = [prototype(PUBLISHED["shelf"]), prototype(PUBLISHED["high_pass"])]


def biquad(p: dict, sr: int) -> list:
    """The bilinear transform at sr, prewarped at f0: s -> (1/K)(1 - 1/z)/(1 + 1/z), K = tan(pi f0 / sr)."""
    k = mp.tan(mp.pi * p["f0"] / sr)
    kq, k2 = k / p["q"], k * k
    d = 1 + kq + k2
    return [
        (p["vh"] + p["vb"] * kq + p["vl"] * k2) / d,
        2 * (p["vl"] * k2 - p["vh"]) / d,
        (p["vh"] - p["vb"] * kq + p["vl"] * k2) / d,
        2 * (k2 - 1) / d,
        (1 - kq + k2) / d,
    ]


def k_weighting(sr: int) -> list[list[float]]:
    return [[f(c) for c in biquad(p, sr)] for p in PROTOTYPES]


def k_weight(x: np.ndarray, sr: int) -> np.ndarray:
    """x of shape (n,) or (n, channels) through both stages, each channel from rest."""
    y = np.asarray(x, dtype=float)
    for b0, b1, b2, a1, a2 in k_weighting(sr):
        y = scipy.signal.lfilter([b0, b1, b2], [1.0, a1, a2], y, axis=0)
    return y


def filter_vectors() -> None:
    claim = "safety:K-WEIGHTING-COEFFICIENTS"
    published = [[f(c) for c in PUBLISHED[s]] for s in ("shelf", "high_pass")]
    cases = [
        case(claim, "k_weighting", [48000], published, "the published BS.1770 coefficients", {"rel": 0, "abs": 1e-13})
    ]
    cases += [
        case(claim, "k_weighting", [sr], k_weighting(sr))
        for sr in (8000, 16000, 22050, 32000, 44100, 48000, 96000, 192000)
    ]
    cases.append(
        case(
            claim,
            "k_weighting",
            [3364],
            k_weighting(3364),
            "the lowest rate: tan is ill-conditioned",
            {"rel": 1e-9, "abs": 0},
        )
    )
    impulse = np.zeros(48)
    impulse[0] = 1.0
    step = np.ones(40)
    stereo = interleave(tones(44100, 0.0008, [(0.5, 3000.0, 0.2)], 6), tones(44100, 0.0008, [(0.4, 15000.0, 1.0)], 6))
    for x, channels, sr, note in (
        (impulse, 1, 48000, "the impulse response at 48 kHz"),
        (step, 1, 8000, "the step response at 8 kHz"),
        (stereo, 2, 44100, "stereo at 44.1 kHz: channels filtered separately"),
    ):
        y = k_weight(x.reshape(-1, channels), sr).reshape(-1)
        cases.append(case("safety:K-WEIGHT-FILTER", "k_weight", [x.tolist(), channels, sr], y.tolist(), note))
    write("filters", "K-weighting coefficients (mpmath, 50 digits) and the filter (scipy.signal.lfilter).", cases)


# Loudness: pyloudnorm's meter with the reference filter, and a transcription ---------------


class Stage:
    """A fixed biquad with pyloudnorm's filter interface."""

    def __init__(self, c: list[float]) -> None:
        self.b, self.a = c[:3], [1.0, *c[3:]]

    def apply_filter(self, data: np.ndarray) -> np.ndarray:
        return scipy.signal.lfilter(self.b, self.a, data)


def meter(sr: int, block: float = 0.4, overlap: float = 0.75) -> pyloudnorm.Meter:
    m = pyloudnorm.Meter(sr, filter_class="custom", block_size=block, overlap=overlap)
    m._filters = {"shelf": Stage(k_weighting(sr)[0]), "high_pass": Stage(k_weighting(sr)[1])}
    return m


def bs1770(frames: np.ndarray, weights: list[float], sr: int, hops: int) -> tuple[list[float], float | None]:
    """BS.1770 from the text: window loudness every 100 ms over `hops` hops, and the gated integrated loudness."""
    y = k_weight(frames, sr)
    n = len(y)
    bound = [k * sr // 10 for k in range(10 * n // sr + 1)]
    powers = []
    for j in range(len(bound) - hops):
        z = np.mean(np.square(y[bound[j] : bound[j + hops]]), axis=0)
        powers.append(float(np.dot(weights, z)))
    with np.errstate(divide="ignore"):
        level = [-0.691 + 10 * np.log10(p) for p in powers]
    above = [p for p, lv in zip(powers, level, strict=True) if lv > -70]
    if not above:
        return level, None
    gate = -0.691 + 10 * np.log10(np.mean(above)) - 10
    gated = [p for p, lv in zip(powers, level, strict=True) if lv > max(gate, -70)]
    return level, -0.691 + 10 * np.log10(np.mean(gated))


def loudness_vectors() -> None:
    claim = "safety:LOUDNESS-MATCHES-PYLOUDNORM"
    a, b, c = (0.3, 997.0, 0.0), (0.1, 60.0, 1.0), (0.05, 1800.0, 2.0)
    loud = tones(4000, 0.8, [(0.3, 997.0, 0.0), (0.2, 150.0, 0.5)])
    middle = tones(4000, 0.6, [(0.006, 440.0, 0.0), (0.004, 1500.0, 1.0)])
    faint = tones(4000, 0.6, [(0.0003, 1000.0, 0.0)])
    five = [tones(4000, 0.4, [(0.1 * (k + 1), 200.0 * (k + 1), k)]) for k in range(5)]
    six = [tones(4000, 0.4, [(0.08 * (k + 1), 150.0 * (k + 1) + 7, 0.5 * k)]) for k in range(6)]
    signals = [
        ("mono, steady", tones(4000, 0.5, [a, b, c]), [1.0], 4000, True),
        (
            "mono: loud, then about 34 dB quieter (below the relative gate), then below the absolute gate",
            np.concatenate([loud, middle, faint]),
            [1.0],
            4000,
            True,
        ),
        (
            "stereo",
            interleave(tones(4000, 0.6, [a, b]), tones(4000, 0.6, [(0.2, 300.0, 0.0), c])),
            [1.0, 1.0],
            4000,
            True,
        ),
        ("5.0 surround: L, R, C, Ls, Rs", interleave(*five), [1.0, 1.0, 1.0, 1.41, 1.41], 4000, True),
        (
            "5.1 with LFE (weight 0); reference: bs1770 only",
            interleave(*six),
            [1.0, 1.0, 1.0, 0.0, 1.41, 1.41],
            4000,
            False,
        ),
        ("mono at 8 kHz, 0.45 s: one block and half a hop", tones(8000, 0.45, [a, c]), [1.0], 8000, True),
    ]
    cases = []
    for note, x, weights, sr, by_pyloudnorm in signals:
        frames = x.reshape(-1, len(weights))
        _, reference = bs1770(frames, weights, sr, 4)
        if by_pyloudnorm:
            got = meter(sr).integrated_loudness(frames if len(weights) > 1 else x)
            assert abs(got - reference) < 1e-11, (note, got, reference)
            reference = got
        cases.append(case(claim, "integrated_loudness", [x.tolist(), weights, sr], float(reference), note))
    bursts = np.concatenate([tones(4000, 0.3, [a]), tones(4000, 0.4, [(0.02, 500.0, 0.0)]), tones(4000, 0.3, [a, b])])
    m = meter(4000)
    m.integrated_loudness(bursts)
    level, _ = bs1770(bursts.reshape(-1, 1), [1.0], 4000, 4)
    assert np.allclose(m.blockwise_loudness, level, rtol=0, atol=1e-11)
    cases.append(
        case(claim, "momentary_loudness", [bursts.tolist(), [1.0], 4000], m.blockwise_loudness, "1 s: 7 windows")
    )
    # Short-term: bs1770 is the reference. pyloudnorm's meter with 3 s blocks and overlap 1 - 1/30 computes its
    # block bounds as int(3.0 * (j * (1 - overlap)) * sr), which falls one sample short (399 for 400 at j = 1).
    long = np.concatenate([tones(4000, 1.6, [a, b]), tones(4000, 1.7, [(0.05, 700.0, 0.0), c])])
    level, _ = bs1770(long.reshape(-1, 1), [1.0], 4000, 30)
    note = "3.3 s: 4 windows; reference: bs1770 only"
    cases.append(case(claim, "short_term_loudness", [long.tolist(), [1.0], 4000], [float(v) for v in level], note))
    write(
        "loudness",
        "Integrated, momentary and short-term loudness: pyloudnorm with the reference K-weighting.",
        cases,
        LU,
    )


# EBU Tech 3341 test signals for tests/ebu3341.rs ---------------------------------------------

EBU_SR, EBU_HZ = 48000, 1000.0
STEREO, FIVE = [1.0, 1.0], [1.0, 1.0, 1.0, 1.41, 1.41]
EBU_CASES = [
    # case, measure, weights, segments [(seconds, dBFS per channel)], expected LUFS; tolerance +-0.1 LU throughout
    (1, "integrated", STEREO, [(20.0, -23.0)], -23.0),
    (1, "momentary", STEREO, [(20.0, -23.0)], -23.0),
    (1, "short_term", STEREO, [(20.0, -23.0)], -23.0),
    (2, "integrated", STEREO, [(20.0, -33.0)], -33.0),
    (2, "momentary", STEREO, [(20.0, -33.0)], -33.0),
    (2, "short_term", STEREO, [(20.0, -33.0)], -33.0),
    (3, "integrated", STEREO, [(10.0, -36.0), (60.0, -23.0), (10.0, -36.0)], -23.0),
    (4, "integrated", STEREO, [(10.0, -72.0), (10.0, -36.0), (60.0, -23.0), (10.0, -36.0), (10.0, -72.0)], -23.0),
    (5, "integrated", STEREO, [(20.0, -26.0), (20.1, -20.0), (20.0, -26.0)], -23.0),
    (6, "integrated", FIVE, [(20.0, [-28.0, -28.0, -24.0, -30.0, -30.0])], -23.0),
    (9, "short_term", STEREO, [(1.34, -20.0), (1.66, -30.0)] * 5, -23.0),
]


def ebu_signal(weights: list[float], segments) -> np.ndarray:
    """A 1 kHz sine, sample i = sin(2 pi f (i + 1) / sr) as signal_core::osc, at each segment's level."""
    n = [round(sec * EBU_SR) for sec, _ in segments]
    tone = np.sin(2 * np.pi * np.cumsum(np.full(sum(n), EBU_HZ)) / EBU_SR)
    gains = [np.broadcast_to(10 ** (np.asarray(dbfs, dtype=float) / 20), len(weights)) for _, dbfs in segments]
    return tone[:, None] * np.repeat(np.array(gains), n, axis=0)


def ebu_reference() -> None:
    out = []
    for number, measure, weights, segments, expected in EBU_CASES:
        entry = {
            "case": number,
            "measure": measure,
            "weights": weights,
            "segments": [[sec, dbfs] for sec, dbfs in segments],
            "expected": expected,
            "tolerance": 0.1,
        }
        if measure == "integrated":
            entry["pyloudnorm"] = float(
                pyloudnorm.Meter(EBU_SR, filter_class="DeMan").integrated_loudness(ebu_signal(weights, segments))
            )
        out.append(entry)
    doc = {
        "description": "EBU Tech 3341 test cases (https://tech.ebu.ch/docs/tech/tech3341.pdf): stereo or 5.0 sines of "
        "1 kHz at 48 kHz, with the published expected loudness and tolerance (+-0.1 LU), and pyloudnorm's integrated "
        "loudness of the same signals (DeMan filters, which reproduce the BS.1770 coefficients at 48 kHz).",
        "sr": EBU_SR,
        "hz": EBU_HZ,
        "cases": out,
    }
    (TESTS / "ebu3341.json").write_text(json.dumps(doc, indent=1) + "\n", encoding="utf-8", newline="\n")


# Invalid input ---------------------------------------------------------------------------------


def error_vectors() -> None:
    claim = "safety:INPUT-VALIDATION"
    tone = tones(4000, 0.001, [(0.3, 997.0, 0.0)]).tolist()  # 4 samples: these checks come before any length check
    cases = [
        error(claim, "amplitude_to_db", [0.0], "silence has no level"),
        error(claim, "amplitude_to_db", [-1.0], "a negative ratio"),
        error(claim, "power_to_db", [0.0], "silence has no level"),
        error(claim, "db_to_power", [4000.0], "the ratio overflows an f64"),
        error(claim, "db_to_amplitude", [7000.0], "the ratio overflows an f64"),
        error(claim, "sum_incoherent_db", [[]], "no levels"),
        error(claim, "energetic_mean_db", [[]], "no levels"),
        error(claim, "leq", [[60.0], [0.0]], "no time at any level"),
        error(claim, "leq", [[60.0, 70.0], [1.0]], "durations and levels differ in length"),
        error(claim, "leq", [[60.0], [-1.0]], "a negative duration"),
        error(claim, "rms_dbfs", [[]], "an empty signal"),
        error(claim, "rms_dbfs", [[0.0, 0.0, 0.0]], "silence"),
        error(claim, "peak_dbfs", [[0.0]], "silence"),
        error(claim, "normalize", [[0.25] * 64, 1, -20.0, -1.0, 1000], "a constant is silent after the high-pass"),
        error(
            claim,
            "normalize",
            [[0.3, -0.2] * 500, 2, -20.0, -1.0, 1000],
            "a constant in each channel, not exact in binary: only rounding residue is left after the high-pass",
        ),
        error(claim, "normalize", [[0.1, 0.2, 0.3], 2, -20.0, -1.0, 1000], "not a whole number of stereo frames"),
        error(claim, "normalize", [[0.1, 0.2], 1, -20.0, -1.0, 0], "sr = 0"),
        error(claim, "allowed_time", [85.0, 85.0, 0.0, 8.0], "an exchange rate of 0"),
        error(claim, "allowed_time", [85.0, 85.0, -3.0, 8.0], "a negative exchange rate"),
        error(claim, "allowed_time", [85.0, 85.0, 3.0, 0.0], "a criterion duration of 0"),
        error(claim, "allowed_time", [-5000.0, 85.0, 3.0, 8.0], "the allowed time overflows an f64"),
        error(claim, "noise_dose", [[90.0], [1.0, 2.0], 85.0, 3.0, 8.0], "hours and levels differ in length"),
        error(claim, "noise_dose", [[90.0], [-1.0], 85.0, 3.0, 8.0], "a negative duration"),
        error(claim, "twa_from_dose", [0.0, 85.0, 3.0], "a dose of 0 has no level"),
        error(claim, "twa_from_dose", [-50.0, 85.0, 3.0], "a negative dose"),
        error(claim, "k_weighting", [3363], "below the lowest rate (the shelf must lie below Nyquist)"),
        error(claim, "k_weighting", [-48000], "a negative sample rate"),
        error(claim, "k_weight", [[0.1, 0.2, 0.3], 2, 48000], "not a whole number of stereo frames"),
        error(claim, "k_weight", [[0.1], 0, 48000], "no channels"),
        error(claim, "k_weight", [[1e308, -1e308] * 4, 1, 48000], "the filter output overflows an f64"),
        error(
            claim,
            "integrated_loudness",
            [tones(4000, 0.39975, [(0.3, 997.0, 0.0)], 2).tolist(), [1.0], 4000],
            "1599 samples: shorter than one 400 ms block",
        ),
        error(
            claim, "integrated_loudness", [[0.0] * 1600, [1.0], 4000], "silence: every block below the absolute gate"
        ),
        error(claim, "integrated_loudness", [[1e200] * 1600, [1.0], 4000], "the K-weighted energy overflows an f64"),
        error(claim, "integrated_loudness", [tone, [], 4000], "no channel weights"),
        error(claim, "integrated_loudness", [tone, [-1.0], 4000], "a negative channel weight"),
        error(claim, "integrated_loudness", [tone[:3], [1.0, 1.0], 4000], "not a whole number of stereo frames"),
        error(claim, "momentary_loudness", [tone, [1.0], 3000], "below the lowest rate"),
        error(claim, "short_term_loudness", [tone, [1.0], 0], "sr = 0"),
    ]
    write("errors", "Invalid input is rejected with an error in every binding (the crate's own input rules).", cases)


if __name__ == "__main__":
    decibel_vectors()
    dose_vectors()
    level_vectors()
    filter_vectors()
    loudness_vectors()
    error_vectors()
    ebu_reference()
    for p in PROTOTYPES:
        print({k: repr(float(v)) for k, v in p.items()})
