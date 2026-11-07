// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2026 Chethan Dandebala Venkatarama

//! Golden-vector tests. The Python and WebAssembly bindings check the same
//! `vectors/*.json` cases, so all three implementations give the same answers.

use audio_safety_math as asm;
use golden_vectors::{Outcome, f64_arg, i64_arg};
use serde_json::{Value, json};
use std::error::Error;

type Res = Result<Value, Box<dyn Error>>;

/// Argument `i` as a whole number of type `T`; out-of-range values (a negative
/// sample rate, say) are rejected, as the Python and JavaScript bindings reject them.
fn int<T: TryFrom<i64>>(args: &[Value], i: usize) -> Result<T, Box<dyn Error>>
where
    T::Error: Error + 'static,
{
    Ok(T::try_from(i64_arg(args, i))?)
}

fn floats(args: &[Value], i: usize) -> Vec<f64> {
    let items = args[i].as_array().expect("an array of numbers");
    items
        .iter()
        .map(|v| v.as_f64().expect("a number"))
        .collect()
}

/// Runs the crate's function `name`; `None` if there is no such function.
fn call(name: &str, a: &[Value]) -> Option<Res> {
    let x = |i| f64_arg(a, i);
    let v = |i| floats(a, i);
    let run = |f: &dyn Fn() -> Res| Some(f());
    match name {
        "amplitude_to_db" => run(&|| Ok(json!(asm::amplitude_to_db(x(0))?))),
        "db_to_amplitude" => run(&|| Ok(json!(asm::db_to_amplitude(x(0))?))),
        "power_to_db" => run(&|| Ok(json!(asm::power_to_db(x(0))?))),
        "db_to_power" => run(&|| Ok(json!(asm::db_to_power(x(0))?))),
        "sum_incoherent_db" => run(&|| Ok(json!(asm::sum_incoherent_db(&v(0))?))),
        "energetic_mean_db" => run(&|| Ok(json!(asm::energetic_mean_db(&v(0))?))),
        "leq" => run(&|| Ok(json!(asm::leq(&v(0), &v(1))?))),
        "rms_dbfs" => run(&|| Ok(json!(asm::rms_dbfs(&v(0))?))),
        "peak_dbfs" => run(&|| Ok(json!(asm::peak_dbfs(&v(0))?))),
        "normalize" => run(&|| {
            Ok(json!(asm::normalize(
                &v(0),
                int(a, 1)?,
                x(2),
                x(3),
                int(a, 4)?
            )?))
        }),
        "allowed_time" => run(&|| Ok(json!(asm::allowed_time(x(0), x(1), x(2), x(3))?))),
        "noise_dose" => run(&|| Ok(json!(asm::noise_dose(&v(0), &v(1), x(2), x(3), x(4))?))),
        "twa_from_dose" => run(&|| Ok(json!(asm::twa_from_dose(x(0), x(1), x(2))?))),
        "k_weighting" => run(&|| Ok(json!(asm::k_weighting(int(a, 0)?)?))),
        "k_weight" => run(&|| Ok(json!(asm::k_weight(&v(0), int(a, 1)?, int(a, 2)?)?))),
        "integrated_loudness" => {
            run(&|| Ok(json!(asm::integrated_loudness(&v(0), &v(1), int(a, 2)?)?)))
        }
        "momentary_loudness" => {
            run(&|| Ok(json!(asm::momentary_loudness(&v(0), &v(1), int(a, 2)?)?)))
        }
        "short_term_loudness" => {
            run(&|| Ok(json!(asm::short_term_loudness(&v(0), &v(1), int(a, 2)?)?)))
        }
        _ => None,
    }
}

/// Maps each vector's `fn` name to this crate's implementation.
fn dispatch(name: &str, args: &[Value]) -> Outcome {
    call(name, args).map_or_else(|| Outcome::unknown(name, args), Outcome::from)
}

#[test]
fn golden_vectors() {
    golden_vectors::check_dir(concat!(env!("CARGO_MANIFEST_DIR"), "/vectors"), dispatch);
}
