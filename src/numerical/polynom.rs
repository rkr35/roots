// Copyright (c) 2017, Mikhail Vorotilov
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

use super::super::find_roots_cubic;
use super::super::find_roots_linear;
use super::super::find_roots_quadratic;
use super::super::find_roots_quartic;
use super::super::FloatType;
use super::super::analytical::roots::Roots;

use super::Convergency;
use super::Interval;
use super::Sample;
use super::SearchError;

#[derive(Debug, PartialEq)]
struct ValueAndDerivative<F>
where
    F: FloatType,
{
    value: Sample<F>,
    derivative: F,
}

trait Polynom<F>
where
    F: FloatType,
{
    fn value(&self, x: &F) -> F;
    fn value_and_derivative(&self, x: &F) -> ValueAndDerivative<F>;
    fn find_root(&self, bracketed_start: &mut Interval<F>, convergency: &mut Convergency<F>) -> Result<F, SearchError>;
}

impl<F> Polynom<F> for [F]
where
    F: FloatType,
{
    fn value(&self, x: &F) -> F {
        let mut result = F::zero();
        let mut xn = F::one();

        // Sum starting with a0
        for a in self.iter().rev() {
            result = result + *a * xn;
            xn = xn * *x;
        }

        // The highest coefficient of the normalized polynom is 1
        result + xn
    }

    fn value_and_derivative(&self, x: &F) -> ValueAndDerivative<F> {
        let mut xn = F::one(); // x^n for SUM(A(n)*x^(n))
        let mut value = F::zero();

        let mut xn1 = F::zero(); // x^n-1 for SUM(n*A(n-1)*x^(n-1))
        let mut derivative = F::zero();
        let mut n = F::zero();

        // Sum starting with a0
        for a in self.iter().rev() {
            value = value + *a * xn;
            derivative = derivative + *a * n * xn1;
            xn1 = xn;
            xn = xn * *x;
            n = n + F::one();
        }

        // The highest coefficient of the normalized polynom is 1
        ValueAndDerivative {
            value: Sample { x: *x, y: value + xn },
            derivative: derivative + n * xn1,
        }
    }

    fn find_root(&self, bracketed_start: &mut Interval<F>, convergency: &mut Convergency<F>) -> Result<F, SearchError> {
        if bracketed_start.is_bracketed() {
            let interval = bracketed_start;
            let mut iter = 0;
            loop {
                if convergency.is_root_found(interval.begin.y) {
                    break Ok(interval.begin.x);
                } else if convergency.is_root_found(interval.end.y) {
                    break Ok(interval.end.x);
                } else if interval.is_converged(convergency) {
                    break Ok(interval.middle());
                } else {
                    let middle = self.value_and_derivative(&interval.middle());
                    let next_sample = if middle.derivative != F::zero() {
                        let newton_raphson = middle.value.x - middle.value.y / middle.derivative;
                        if newton_raphson >= interval.begin.x && newton_raphson <= interval.end.x {
                            let newton_raphson_value = self.value(&newton_raphson);
                            if newton_raphson_value.abs() < middle.value.y.abs() {
                                Sample {
                                    x: newton_raphson,
                                    y: newton_raphson_value,
                                }
                            } else {
                                middle.value
                            }
                        } else {
                            middle.value
                        }
                    } else {
                        middle.value
                    };
                    if interval.begin.is_bracketed_with(&next_sample) {
                        interval.end = Sample {
                            x: next_sample.x,
                            y: next_sample.y,
                        };
                    } else {
                        interval.begin = Sample {
                            x: next_sample.x,
                            y: next_sample.y,
                        };
                    }
                }
                iter = iter + 1;
                if convergency.is_iteration_limit_reached(iter) {
                    break Err(SearchError::NoConvergency);
                }
            }
        } else {
            Err(SearchError::NoBracketing)
        }
    }


}


/// Find all roots of the normalized polynomial
/// x^n + a[0]*x^(n-1) + a[1]*x^(n-2) + ... + a[n-1] = 0
/// using the Sturm's theorem recursively.
///
/// # Examples
///
/// ```
/// use roots::find_roots_sturm;
///
/// let polynom = &[1f64,1f64,1f64,1f64,1f64,1f64];
///
/// let roots_or_errors = find_roots_sturm(polynom, &mut 1e-6);
/// // Returns vector of roots or search errors;
///
/// let roots: Vec<_> = find_roots_sturm(polynom, &mut 1e-8f64)
///             .iter()
///             .filter_map(|s| match s {
///                 &Ok(ref x) => Some(*x),
///                 &Err(_) => None,
///             })
///             .collect();
/// // Returns vector of roots filterin out all search errors;
/// ```
pub fn find_roots_sturm<F>(a: &[F]) -> Option<Roots<F>>
where
    F: FloatType,
{
    Some(match a.len() {
        0 => Roots::No([]),
        1 => find_roots_linear(F::one(), a[0]),
        2 => find_roots_quadratic(F::one(), a[0], a[1]),
        3 => find_roots_cubic(F::one(), a[0], a[1], a[2]),
        4 => find_roots_quartic(F::one(), a[0], a[1], a[2], a[3]),
        _ => {
            return None;
        },
    })
}

#[cfg(test)]
mod test {
    use super::super::*;
    use super::*;

    #[test]
    fn test_find_roots_sturm() {
        let polynom = &[-2f64, 1f64];
        let roots = find_roots_sturm(polynom, &mut 1e-6f64);
        assert_eq!(roots, [Ok(1f64)]);
    }

    #[test]
    fn test_polynom_value() {
        let polynom = [1f64, -2f64, 1f64];
        assert_eq!(1f64, polynom.value(&0f64));
        assert_eq!(1f64, polynom.value(&1f64));
        assert_eq!(3f64, polynom.value(&-1f64));
    }

    #[test]
    fn test_polynom_value_and_derivative() {
        let polynom = [1f64, -2f64, 1f64];
        assert_eq!(
            ValueAndDerivative {
                value: Sample { x: 0f64, y: 1f64 },
                derivative: -2f64
            },
            polynom.value_and_derivative(&0f64)
        );
        assert_eq!(
            ValueAndDerivative {
                value: Sample { x: 1f64, y: 1f64 },
                derivative: 3f64
            },
            polynom.value_and_derivative(&1f64)
        );
        assert_eq!(
            ValueAndDerivative {
                value: Sample { x: -1f64, y: 3f64 },
                derivative: -1f64
            },
            polynom.value_and_derivative(&-1f64)
        );
    }

    #[test]
    fn test_derivative_polynom_3() {
        // x^3 + 1*x^2 - 2*x^1 + 1*x^0 => 3*x^2 + 2*x^1 - 2*x^0 => x^2 + (2/3)*x^1 - (2/3)*x^0
        let polynom = [1f64, -2f64, 1f64];
        let derivative = polynom.derivative_polynom();
        assert_float_array_eq!(1e-15, derivative, [2f64 / 3f64, -2f64 / 3f64]);
    }

    #[test]
    fn test_derivative_polynom_5() {
        // x^5 - 2*x^4 - 3*x^3 + 4*x^2 + 0*x^1 + 0*x^0 => 5*x^4 - 8*x^3 - 9*x^2 + 8*x^1 + 0*x^0 => x^4 - (8/5)*x^3 - (9/5)*x^2 + (8/5)*x^1 + 0*x^0
        let polynom = [-2f64, -3f64, 4f64, 0f64, 0f64];
        let derivative = polynom.derivative_polynom();
        assert_float_array_eq!(1e-15, derivative, [-8f64 / 5f64, -9f64 / 5f64, 8f64 / 5f64, 0f64]);
    }

    #[test]
    fn find_roots_sturm_7() {
        // x^7+4.0*x^6-4.0*x^4+2.0*x^3+1.0*x^2+6.0*x^1-3.0*x^0 => {-3.6547, -1.67175, 0.455904}
        let polynom = [4f64, 0f64, -4f64, 2f64, 1f64, 6f64, -3f64];
        let roots: Vec<_> = find_roots_sturm(&polynom, &mut 1e-8f64)
            .iter()
            .filter_map(|s| match s {
                &Ok(ref x) => Some(*x),
                &Err(_) => None,
            })
            .collect();
        assert_float_array_eq!(1e-5, roots, [-3.6547f64, -1.67175f64, 0.455904f64]);
    }
}
