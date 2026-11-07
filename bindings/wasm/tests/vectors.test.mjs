// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
// The shared golden vectors, checked against the WebAssembly build.

import assert from "node:assert/strict";
import { createRequire } from "node:module";
import { test } from "node:test";
import { checkDir } from "golden-vectors";

const impl = createRequire(import.meta.url)("../pkg/node/audio_safety_math_wasm.js");

test("golden vectors", () => {
  assert.ok(checkDir(new URL("../../../vectors/", import.meta.url), impl) > 0);
});

test("dose functions default to the NIOSH rule", () => {
  assert.equal(impl.allowed_time(88), 4);
  assert.equal(impl.allowed_time(95, 90, 5), 4);
  assert.equal(impl.noise_dose([88, 85], [2, 4]), 100);
  assert.equal(impl.twa_from_dose(200), 88);
});

test("signals are typed arrays, silence is -Infinity in the meters, errors throw", () => {
  const tone = Float64Array.from({ length: 8000 }, (_, i) => 0.1 * Math.sin((2 * Math.PI * 997 * i) / 8000));
  assert.ok(impl.k_weight(tone, 1, 8000) instanceof Float64Array);
  const [shelf, highPass] = impl.k_weighting(48000);
  assert.ok(shelf instanceof Float64Array && highPass.length === 5);
  const gap = new Float64Array(8000);
  gap.set(tone.subarray(0, 4000), 4000);
  const m = impl.momentary_loudness(gap, [1], 8000);
  assert.equal(m[0], -Infinity);
  assert.ok(Number.isFinite(m[m.length - 1]));
  assert.throws(() => impl.k_weighting(48000.5), /sr must be a whole number/);
  assert.throws(() => impl.k_weighting(-1), /sr must be a whole number/);
  assert.throws(() => impl.amplitude_to_db(0), /has no level in decibels/);
  assert.throws(() => impl.integrated_loudness(new Float64Array(4000), [1], 8000), /absolute gate/);
});
