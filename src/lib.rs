// Copyright (c) 2015, Mikhail Vorotilov
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
//
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// #![no_std]
//#![crate_id = "roots"]
#![crate_type = "lib"]

//! A set of functions to find real roots of numerical equations.
//!
//! This crate contains various algorithms for numerical and analytical solving
//! of 1-variable equations like f(x)=0. Only real roots are calculated.
//! Multiple (double etc.) roots are considered as one root.
//!
//! # Use
//!
//! Functions find_root_* try to find a root of any given closure function by
//! iterative approximations. Conditions for success/failure can be customized
//! by implementing the Convergency trait.
//! Functions find_roots_* return all roots of several simple equations at once.

#[cfg(test)]
macro_rules! assert_float_eq(
    ($precision:expr, $given:expr , $expected:expr) => ({
      match (&($precision), &($given), &($expected)) {
          (precision_val, given_val, expected_val) => {
            let diff = given_val-expected_val;
            if diff.abs() > precision_val.abs() {
              panic!("floats are not the same: (`{}`: `{:.15e}`, expected: `{:.15e}`, precision: `{:.15e}`, delta: `{:.15e}`)", stringify!($given), *given_val, *expected_val, *precision_val, diff )
            }
          }
      }
    })
);

#[cfg(test)]
macro_rules! assert_float_array_eq(
    ($precision:expr, $given:expr , $expected:expr) => ({

      let expected_len = $expected.len();
      let mut seen = 0;
      $given
        .zip($expected.iter())
        .for_each(|(g, &e)| {
          seen += 1;
          assert_float_eq!($precision, g, e);
        });

      assert_eq!(seen, expected_len);

    })
);

pub mod analytical;
pub mod float;

pub use self::float::FloatType;

pub use self::analytical::biquadratic::find_roots_biquadratic;
pub use self::analytical::cubic::find_roots_cubic;
pub use self::analytical::cubic_depressed::find_roots_cubic_depressed;
pub use self::analytical::cubic_normalized::find_roots_cubic_normalized;
pub use self::analytical::linear::find_roots_linear;
pub use self::analytical::quadratic::find_roots_quadratic;
pub use self::analytical::quartic::find_roots_quartic;
pub use self::analytical::quartic_depressed::find_roots_quartic_depressed;
pub use self::analytical::roots::Roots;
