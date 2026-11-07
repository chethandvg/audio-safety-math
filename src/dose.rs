// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Noise dose by the exchange-rate rule: allowed time, daily dose and the
//! equivalent eight-hour level.

use crate::{Result, all_finite, all_non_negative, finite, fits, positive, same_length};

/// The NIOSH recommended exposure limit: 85 dBA as an 8-hour time-weighted
/// average (NIOSH 1998, publication 98-126, Chapter 1, Section 1.1,
/// doi:10.26616/NIOSHPUB98126; `safety:NIOSH-REL`, level R).
pub const NIOSH_CRITERION_DB: f64 = 85.0;

/// The NIOSH exchange rate: every 3 dB more halves the allowed time (same
/// source as [`NIOSH_CRITERION_DB`]).
pub const NIOSH_EXCHANGE_DB: f64 = 3.0;

/// The criterion duration, 8 hours (a working day).
pub const CRITERION_HOURS: f64 = 8.0;

/// The criterion level and exchange rate, checked.
fn rule(criterion_db: f64, exchange_db: f64) -> Result<(f64, f64)> {
    Ok((
        finite("criterion_db", criterion_db)?,
        positive("exchange_db", exchange_db)?,
    ))
}

/// `(L − L₀)/q`: how many times the allowed time halves between the criterion
/// level and `level`, so that an hour at `level` is worth `2^doublings` hours
/// at the criterion level.
fn doublings(level: f64, criterion_db: f64, exchange_db: f64) -> f64 {
    (level - criterion_db) / exchange_db
}

/// The time allowed at `level_dba` before the daily dose reaches 100 %:
/// `T = criterion_hours / 2^((L − criterion_db) / exchange_db)`, in the unit of
/// `criterion_hours`. With the NIOSH parameters (85 dBA, 3 dB, 8 h) 88 dBA is
/// allowed 4 hours and 100 dBA 15 minutes; the OSHA rule of 90 dBA with a 5 dB
/// exchange rate is `allowed_time(L, 90.0, 5.0, 8.0)`. Raising the level by
/// one exchange rate halves the time (proved as `safety:EXCHANGE-RATE`). A time
/// too small for an `f64` rounds to 0. When `2^((L − criterion_db) / exchange_db)`
/// is outside the normal range of an `f64` (about 1022 exchange rates or more
/// from the criterion level), the time is computed instead as the one power of
/// two `2^(log₂ criterion_hours − (L − criterion_db) / exchange_db)`, so a time
/// that fits in an `f64` is not lost on the way: with the NIOSH parameters
/// 3160 dBA is allowed 2⁻¹⁰²² hours.
pub fn allowed_time(
    level_dba: f64,
    criterion_db: f64,
    exchange_db: f64,
    criterion_hours: f64,
) -> Result<f64> {
    let level = finite("level_dba", level_dba)?;
    let (l0, q) = rule(criterion_db, exchange_db)?;
    let tc = positive("criterion_hours", criterion_hours)?;
    let d = doublings(level, l0, q);
    let p = d.exp2();
    fits(if p.is_normal() {
        tc / p
    } else {
        (tc.log2() - d).exp2()
    })
}

/// The daily noise dose in percent: `D = 100 · Σ Cᵢ / Tᵢ`, where segment `i`
/// lasts `hours[i]` (`Cᵢ`, in the unit of `criterion_hours`) at
/// `levels_dba[i]` and `Tᵢ` is its [`allowed_time`]. Each share is computed
/// as `Cᵢ / Tᵢ = 2^eᵢ` with `eᵢ = log₂ Cᵢ − log₂ criterion_hours +
/// (Lᵢ − criterion_db) / exchange_db`, and the shares are summed around the
/// largest exponent, so no division by a tiny allowed time is needed and a
/// short time at a huge level does not overflow on the way: only a dose too
/// large for an `f64` is an error. The largest share, `2^max eᵢ`, is applied
/// last, so a dose whose largest share is below the smallest normal `f64`
/// (2⁻¹⁰²²) has only that share's subnormal precision, and is 0 when the share
/// rounds to 0 (about 2⁻¹⁰⁷⁵ or less), although the dose, at least 100 times the
/// share, may fit in an `f64`. 100 % is the limit; no segments give 0 %, and a
/// segment of 0 hours counts for nothing, however loud it is.
pub fn noise_dose(
    levels_dba: &[f64],
    hours: &[f64],
    criterion_db: f64,
    exchange_db: f64,
    criterion_hours: f64,
) -> Result<f64> {
    all_finite("levels_dba", levels_dba)?;
    same_length("hours", hours, levels_dba.len())?;
    all_non_negative("hours", hours)?;
    let (l0, q) = rule(criterion_db, exchange_db)?;
    let tc = positive("criterion_hours", criterion_hours)?;
    let exponents = || {
        levels_dba
            .iter()
            .zip(hours)
            .filter(|&(_, &c)| c > 0.0)
            .map(|(&l, &c)| (c.log2() - tc.log2()) + doublings(l, l0, q))
    };
    let top = exponents().fold(f64::NEG_INFINITY, f64::max);
    if top.is_infinite() {
        // No time at any level (0 %), or a dose beyond any f64 (an error).
        return fits(100.0 * top.exp2());
    }
    let shares: f64 = exponents().map(|e| (e - top).exp2()).sum();
    fits(100.0 * shares * top.exp2())
}

/// The time-weighted average level of a daily dose: the steady level that
/// would give `dose_percent` over the criterion duration,
/// `TWA = criterion_db + exchange_db · log₂(D / 100)`. It is the exact inverse
/// of the dose rule (proved as `safety:TWA-INVERTS-DOSE`): 100 % is the
/// criterion level, 200 % one exchange rate above it. The form
/// `criterion_db + 10 · log₁₀(D / 100)` treats a 3 dB exchange rate as exactly
/// 10 · log₁₀ 2 dB; it differs from this one by 0.0103 dB per doubling of the
/// dose (95.0 against 94.966 dBA at 1000 % with the NIOSH parameters).
pub fn twa_from_dose(dose_percent: f64, criterion_db: f64, exchange_db: f64) -> Result<f64> {
    let d = positive("dose_percent", dose_percent)?;
    let (l0, q) = rule(criterion_db, exchange_db)?;
    fits(l0 + q * (d / 100.0).log2())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn niosh_and_osha_examples() {
        let niosh = |l| allowed_time(l, 85.0, 3.0, 8.0).unwrap();
        assert_eq!([niosh(85.0), niosh(88.0), niosh(91.0)], [8.0, 4.0, 2.0]);
        assert_eq!(niosh(100.0) * 60.0, 15.0);
        // 2^(3075/3) overflows an f64, 2^(-3285/3) vanishes; the times fit.
        assert_eq!(niosh(3160.0), f64::MIN_POSITIVE); // 2^-1022 hours
        let far = allowed_time(3385.0, 85.0, 3.0, 1e300).unwrap();
        assert!((far / 7.362151829022863e-32 - 1.0).abs() < 1e-12, "{far}");
        let near = allowed_time(-3200.0, 85.0, 3.0, 1e-300).unwrap();
        assert!((near / 4.244682903279331e29 - 1.0).abs() < 1e-12, "{near}");
        assert_eq!(allowed_time(95.0, 90.0, 5.0, 8.0).unwrap(), 4.0);
        assert_eq!(
            noise_dose(&[88.0, 85.0], &[2.0, 4.0], 85.0, 3.0, 8.0).unwrap(),
            100.0
        );
        assert_eq!(noise_dose(&[], &[], 85.0, 3.0, 8.0).unwrap(), 0.0);
        // No time at 5000 dBA is no dose, although 2^(4915/3) overflows an f64.
        assert_eq!(
            noise_dose(&[5000.0, 85.0], &[0.0, 8.0], 85.0, 3.0, 8.0).unwrap(),
            100.0
        );
        // 10^-10 hours at 3160 dBA: 2^(3075/3) overflows an f64, the dose does not.
        let short = noise_dose(&[3160.0], &[1e-10], 85.0, 3.0, 8.0).unwrap();
        let expected = 100.0 / 8.0 * 1e-10 * 2f64.powi(1000) * 2f64.powi(25);
        assert!((short / expected - 1.0).abs() < 1e-12, "{short}");
        assert_eq!(
            noise_dose(&[5000.0], &[8.0], 85.0, 3.0, 8.0)
                .unwrap_err()
                .name(),
            "result"
        );
        assert_eq!(twa_from_dose(200.0, 85.0, 3.0).unwrap(), 88.0);
        assert!((twa_from_dose(1000.0, 85.0, 3.0).unwrap() - 94.966).abs() < 1e-3);
    }

    #[test]
    fn invalid_rules_are_errors() {
        assert_eq!(
            allowed_time(85.0, 85.0, 0.0, 8.0).unwrap_err().name(),
            "exchange_db"
        );
        assert_eq!(
            allowed_time(85.0, 85.0, 3.0, -8.0).unwrap_err().name(),
            "criterion_hours"
        );
        assert_eq!(
            allowed_time(-5000.0, 85.0, 3.0, 8.0).unwrap_err().name(),
            "result"
        );
        assert_eq!(allowed_time(5000.0, 85.0, 3.0, 8.0).unwrap(), 0.0);
        assert!(noise_dose(&[90.0], &[1.0, 2.0], 85.0, 3.0, 8.0).is_err());
        assert!(noise_dose(&[90.0], &[-1.0], 85.0, 3.0, 8.0).is_err());
        assert!(noise_dose(&[f64::NAN], &[1.0], 85.0, 3.0, 8.0).is_err());
        assert!(twa_from_dose(0.0, 85.0, 3.0).is_err());
    }
}
