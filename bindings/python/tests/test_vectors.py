# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
"""The shared golden vectors, checked against the Python bindings."""

import math
from pathlib import Path

import audio_safety_math as asm
import numpy as np
import pytest
from evidence_kernel.vectors import check_dir

VECTORS = Path(__file__).resolve().parents[3] / "vectors"


def plain(value):
    """numpy arrays and tuples as lists, for the vector comparison."""
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, tuple):
        return [plain(v) for v in value]
    return value


def by_name(module):
    """Every public function of `module`, under its own name, returning plain values."""

    def wrap(fn):
        return lambda *args: plain(fn(*args))

    return {name: wrap(fn) for name, fn in vars(module).items() if callable(fn) and not name.startswith("_")}


def test_golden_vectors():
    assert check_dir(VECTORS, by_name(asm)) > 0


def test_dose_functions_default_to_niosh():
    assert asm.allowed_time(88.0) == asm.allowed_time(88.0, 85.0, 3.0, 8.0) == 4.0
    assert asm.allowed_time(95.0, criterion_db=90.0, exchange_db=5.0) == 4.0
    assert asm.noise_dose([88.0, 85.0], [2.0, 4.0]) == 100.0
    assert asm.twa_from_dose(200.0) == 88.0
    assert (asm.NIOSH_CRITERION_DB, asm.NIOSH_EXCHANGE_DB, asm.CRITERION_HOURS) == (85.0, 3.0, 8.0)


def test_signals_are_numpy_arrays_and_meters_show_silence_as_minus_infinity():
    tone = 0.1 * np.sin(2 * np.pi * 997.0 * np.arange(8000) / 8000)
    y = asm.k_weight(tone, 1, 8000)
    assert isinstance(y, np.ndarray) and y.dtype == np.float64 and y.shape == (8000,)
    assert np.array_equal(y, asm.k_weight(tone.tolist(), 1, 8000))
    gap = np.concatenate([np.zeros(4000), tone[:4000]])
    m = asm.momentary_loudness(gap, [1.0], 8000)
    assert m[0] == -math.inf and math.isfinite(m[-1])
    assert asm.integrated_loudness(tone[::2], np.array([1.0]), 4000) == asm.integrated_loudness(
        tone[::2].copy(), [1.0], 4000
    )


def test_errors():
    with pytest.raises(ValueError, match="has no level in decibels"):
        asm.amplitude_to_db(0.0)
    with pytest.raises(ValueError, match="not silent"):
        asm.rms_dbfs([0.0, 0.0])
    with pytest.raises(ValueError, match="absolute gate"):
        asm.integrated_loudness(np.zeros(4000), [1.0], 8000)
    with pytest.raises(ValueError, match="K_WEIGHTING_MIN_SR"):
        asm.k_weighting(3000)
    with pytest.raises(OverflowError):
        asm.k_weighting(-48000)
    with pytest.raises(TypeError):
        asm.k_weighting(48000.5)
    with pytest.raises(TypeError, match="missing"):
        asm.normalize([0.1, 0.2], -20.0)  # lifesound.py's order: the channel count is required
    with pytest.raises(ValueError):
        asm.leq([60.0, float("nan")], [1.0, 1.0])
