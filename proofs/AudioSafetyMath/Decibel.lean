/-
SPDX-License-Identifier: Apache-2.0
SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama
-/
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.SpecialFunctions.Log.Base

/-!
# Decibels

Definitions mirror the Rust functions: `amplitudeDb` is `amplitude_to_db`
(20·log₁₀), `powerDb` is `power_to_db` (10·log₁₀), `dbToPower` is `db_to_power`
and `sumDb` is `sum_incoherent_db`. Everything here is exact real arithmetic.

* `amplitudeDb_mul`: gains in a chain multiply, so their decibels add
  (`safety:DB-PRODUCT`).
* `sumDb_bounds`: the incoherent sum of `n ≥ 1` levels lies between the loudest
  level and that level plus `10·log₁₀ n` (`safety:INCOHERENT-SUM-BOUNDS`).
* `three_db_about_doubling`: doubling the power adds between 3.0102 and
  3.0103 dB (`safety:THREE-DB-IS-ABOUT-DOUBLING`), which is why a 3 dB
  exchange rate approximates equal energy.
-/

namespace AudioSafetyMath

open Real Finset

/-- A power ratio in decibels: `10 · log₁₀ r` (`power_to_db`). -/
noncomputable def powerDb (r : ℝ) : ℝ := 10 * logb 10 r

/-- An amplitude ratio in decibels: `20 · log₁₀ a` (`amplitude_to_db`). -/
noncomputable def amplitudeDb (a : ℝ) : ℝ := 20 * logb 10 a

/-- The power ratio of a level of `L` dB: `10^(L/10)` (`db_to_power`). -/
noncomputable def dbToPower (L : ℝ) : ℝ := (10 : ℝ) ^ (L / 10)

/-- The incoherent sum of the levels `L i` over `s`: `10 · log₁₀ Σ 10^(Lᵢ/10)`
(`sum_incoherent_db`). -/
noncomputable def sumDb {ι : Type*} (s : Finset ι) (L : ι → ℝ) : ℝ :=
  powerDb (∑ i ∈ s, dbToPower (L i))

/-- **safety:DB-PRODUCT.** For amplitude ratios `a, b > 0`,
`20·log₁₀(a·b) = 20·log₁₀ a + 20·log₁₀ b`: gains in a chain add in decibels. -/
theorem amplitudeDb_mul {a b : ℝ} (ha : 0 < a) (hb : 0 < b) :
    amplitudeDb (a * b) = amplitudeDb a + amplitudeDb b := by
  unfold amplitudeDb
  rw [logb_mul ha.ne' hb.ne']
  ring

/-- Decibels and power ratios are inverse: `powerDb (10^(L/10)) = L`. -/
theorem powerDb_dbToPower (L : ℝ) : powerDb (dbToPower L) = L := by
  unfold powerDb dbToPower
  rw [logb_rpow (by norm_num) (by norm_num)]
  ring

theorem dbToPower_pos (L : ℝ) : 0 < dbToPower L := rpow_pos_of_pos (by norm_num) _

/-- `powerDb` is monotone on positive ratios. -/
theorem powerDb_le_powerDb {x y : ℝ} (hx : 0 < x) (hxy : x ≤ y) : powerDb x ≤ powerDb y := by
  unfold powerDb
  have := logb_le_logb_of_le (b := 10) (by norm_num) hx hxy
  linarith

/-- **safety:INCOHERENT-SUM-BOUNDS.** For a nonempty set of `n` levels with
maximum `M`, `M ≤ 10·log₁₀ Σ 10^(Lᵢ/10) ≤ M + 10·log₁₀ n`. -/
theorem sumDb_bounds {ι : Type*} (s : Finset ι) (hs : s.Nonempty) (L : ι → ℝ) :
    s.sup' hs L ≤ sumDb s L ∧ sumDb s L ≤ s.sup' hs L + powerDb s.card := by
  set M := s.sup' hs L with hM
  obtain ⟨i, hi, hiM⟩ := exists_mem_eq_sup' hs L
  have hsum : 0 < ∑ j ∈ s, dbToPower (L j) := sum_pos (fun j _ => dbToPower_pos _) hs
  constructor
  · calc M = powerDb (dbToPower M) := (powerDb_dbToPower M).symm
      _ ≤ sumDb s L := by
        apply powerDb_le_powerDb (dbToPower_pos M)
        rw [hM, hiM]
        exact single_le_sum (fun j _ => (dbToPower_pos (L j)).le) hi
  · have hle : ∑ j ∈ s, dbToPower (L j) ≤ s.card * dbToPower M := by
      have h := sum_le_card_nsmul s (fun j => dbToPower (L j)) (dbToPower M) (fun j hj => by
        have hj' : L j ≤ M := le_sup' L hj
        exact rpow_le_rpow_of_exponent_le (by norm_num) (by linarith))
      simpa [nsmul_eq_mul] using h
    have hcard : (0 : ℝ) < s.card := by exact_mod_cast hs.card_pos
    calc sumDb s L ≤ powerDb (s.card * dbToPower M) := powerDb_le_powerDb hsum hle
      _ = M + powerDb s.card := by
        have h := powerDb_dbToPower M
        unfold powerDb at h ⊢
        rw [logb_mul hcard.ne' (dbToPower_pos M).ne']
        linarith

/-- **safety:THREE-DB-IS-ABOUT-DOUBLING.** `3.0102 < 10·log₁₀ 2 < 3.0103`. The
proof compares whole numbers: `10^30102 < 2^100000 < 10^30103`, checked by the
kernel's arbitrary-precision arithmetic. -/
theorem three_db_about_doubling : 3.0102 < powerDb 2 ∧ powerDb 2 < 3.0103 := by
  have h10 : 0 < log 10 := log_pos (by norm_num)
  have lo : (30102 : ℝ) * log 10 < 100000 * log 2 := by
    have h : (10 : ℕ) ^ 30102 < 2 ^ 100000 := by decide +kernel
    have hr := (Nat.cast_lt (α := ℝ)).mpr h
    rw [Nat.cast_pow, Nat.cast_pow, Nat.cast_ofNat, Nat.cast_ofNat] at hr
    have h' := log_lt_log (by positivity) hr
    rw [log_pow, log_pow] at h'
    exact_mod_cast h'
  have hi : (100000 : ℝ) * log 2 < 30103 * log 10 := by
    have h : (2 : ℕ) ^ 100000 < 10 ^ 30103 := by decide +kernel
    have hr := (Nat.cast_lt (α := ℝ)).mpr h
    rw [Nat.cast_pow, Nat.cast_pow, Nat.cast_ofNat, Nat.cast_ofNat] at hr
    have h' := log_lt_log (by positivity) hr
    rw [log_pow, log_pow] at h'
    exact_mod_cast h'
  unfold powerDb
  rw [logb, mul_div_assoc', lt_div_iff₀ h10, div_lt_iff₀ h10]
  constructor <;> linarith

end AudioSafetyMath
