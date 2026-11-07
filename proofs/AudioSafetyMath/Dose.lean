/-
SPDX-License-Identifier: Apache-2.0
SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
-/
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Log.Base

/-!
# The exchange-rate dose rule

`allowedTime L₀ q T_c L = T_c / 2^((L − L₀)/q)` is the Rust `allowed_time`
(criterion level `L₀`, exchange rate `q`, criterion duration `T_c`), `dose` is a
time divided by its allowed time (the Rust `noise_dose` is 100 times the sum of
these), and `twa` is `twa_from_dose`. Everything here is exact real arithmetic.

* `allowedTime_add_exchange`, `dose_half_time_add_exchange`: one exchange rate
  more halves the allowed time, so `t` hours at `L` and `t/2` hours at `L + q`
  are the same dose (`safety:EXCHANGE-RATE`).
* `twa_dose`: the TWA of a full criterion day at a steady level is that level
  (`safety:TWA-INVERTS-DOSE`).
-/

namespace AudioSafetyMath

open Real

/-- The time allowed at level `L`: `T_c / 2^((L − L₀)/q)` (`allowed_time`). -/
noncomputable def allowedTime (L₀ q Tc L : ℝ) : ℝ := Tc / (2 : ℝ) ^ ((L - L₀) / q)

/-- The dose of `t` hours at level `L`, as a fraction of the limit. -/
noncomputable def dose (L₀ q Tc t L : ℝ) : ℝ := t / allowedTime L₀ q Tc L

/-- The time-weighted average level of a dose `D` in percent:
`L₀ + q · log₂(D / 100)` (`twa_from_dose`). -/
noncomputable def twa (L₀ q D : ℝ) : ℝ := L₀ + q * logb 2 (D / 100)

/-- `k` exchange rates more divide the allowed time by `2^k`, for any real `k`. -/
theorem allowedTime_add_mul (L₀ q Tc L k : ℝ) (hq : q ≠ 0) :
    allowedTime L₀ q Tc (L + k * q) = allowedTime L₀ q Tc L / 2 ^ k := by
  unfold allowedTime
  have h : (L + k * q - L₀) / q = (L - L₀) / q + k := by
    field_simp
    ring
  rw [h, rpow_add (by norm_num : (0 : ℝ) < 2)]
  ring

/-- **safety:EXCHANGE-RATE.** One exchange rate more halves the allowed time:
`allowedTime (L + q) = allowedTime L / 2`. -/
theorem allowedTime_add_exchange (L₀ q Tc L : ℝ) (hq : q ≠ 0) :
    allowedTime L₀ q Tc (L + q) = allowedTime L₀ q Tc L / 2 := by
  have h := allowedTime_add_mul L₀ q Tc L 1 hq
  rwa [one_mul, rpow_one] at h

/-- **safety:EXCHANGE-RATE.** A dose of `t` hours at `L` equals a dose of
`t / 2` hours at `L + q`. -/
theorem dose_half_time_add_exchange (L₀ q Tc t L : ℝ) (hq : q ≠ 0) :
    dose L₀ q Tc (t / 2) (L + q) = dose L₀ q Tc t L := by
  unfold dose
  rw [allowedTime_add_exchange L₀ q Tc L hq]
  ring

/-- **safety:TWA-INVERTS-DOSE.** A whole criterion duration at a steady level
`L` is a dose of `100 · 2^((L − L₀)/q)` percent, whose TWA is `L` again. -/
theorem twa_dose (L₀ q Tc L : ℝ) (hq : q ≠ 0) (hTc : Tc ≠ 0) :
    twa L₀ q (100 * dose L₀ q Tc Tc L) = L := by
  have h2 : (2 : ℝ) ^ ((L - L₀) / q) ≠ 0 := (rpow_pos_of_pos (by norm_num) _).ne'
  have h : 100 * dose L₀ q Tc Tc L / 100 = (2 : ℝ) ^ ((L - L₀) / q) := by
    unfold dose allowedTime
    field_simp
  unfold twa
  rw [h, logb_rpow (by norm_num) (by norm_num)]
  field_simp
  ring

end AudioSafetyMath
