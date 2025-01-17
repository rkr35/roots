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

use core::fmt::Debug;
use core::ops::Add;
use core::ops::Div;
use core::ops::Mul;
use core::ops::Neg;
use core::ops::Sub;

/// Generic type that lists functions and constants needed in calculations.
/// Default implementations for f32 and f64 are provided.
pub trait FloatType:
    Sized
    + Copy
    + Debug
    + From<i16>
    + PartialEq
    + PartialOrd
    + Neg<Output = Self>
    + Add<Output = Self>
    + Sub<Output = Self>
    + Mul<Output = Self>
    + Div<Output = Self>
    + Default
{
    fn zero() -> Self;
    fn one() -> Self;
    fn two() -> Self;
    fn three() -> Self;
    fn pi() -> Self;
    fn one_third() -> Self;
    fn four() -> Self;
    fn five() -> Self {
        Self::two() + Self::three()
    }
    fn nine() -> Self {
        Self::three() * Self::three()
    }
    fn twenty_seven() -> Self {
        Self::nine() * Self::three()
    }
    fn two_third_pi() -> Self;
    fn sqrt(self) -> Self;
    /// The cubic root function is pow(x, 1/3) accepting negative arguments
    fn cbrt(self) -> Self {
        if self < Self::zero() {
            -(-self).powf(Self::one_third())
        } else {
            self.powf(Self::one_third())
        }
    }
    fn acos(self) -> Self;
    fn cos(self) -> Self;
    fn abs(self) -> Self;
    fn powf(self, n: Self) -> Self;
}

impl FloatType for f32 {
    
    fn zero() -> Self {
        0f32
    }
    
    fn one_third() -> Self {
        1f32 / 3f32
    }
    
    fn one() -> Self {
        1f32
    }
    
    fn two() -> Self {
        2f32
    }
    
    fn three() -> Self {
        3f32
    }
    
    fn four() -> Self {
        4f32
    }
    
    fn two_third_pi() -> Self {
        2f32 * core::f32::consts::FRAC_PI_3
    }
    
    fn pi() -> Self {
        core::f32::consts::PI
    }
    fn sqrt(self) -> Self {
        self.sqrt()
    }
    fn acos(self) -> Self {
        self.acos()
    }
    fn cos(self) -> Self {
        self.cos()
    }
    fn abs(self) -> Self {
        self.abs()
    }
    fn powf(self, n: Self) -> Self {
        self.powf(n)
    }
}

impl FloatType for f64 {
    
    fn zero() -> Self {
        0f64
    }
    
    fn one_third() -> Self {
        1f64 / 3f64
    }
    
    fn one() -> Self {
        1f64
    }
    
    fn two() -> Self {
        2f64
    }
    
    fn three() -> Self {
        3f64
    }
    
    fn four() -> Self {
        4f64
    }
    
    fn two_third_pi() -> Self {
        2f64 * core::f64::consts::FRAC_PI_3
    }
    
    fn pi() -> Self {
        core::f64::consts::PI
    }
    fn sqrt(self) -> Self {
        self.sqrt()
    }
    fn acos(self) -> Self {
        self.acos()
    }
    fn cos(self) -> Self {
        self.cos()
    }
    fn abs(self) -> Self {
        self.abs()
    }
    fn powf(self, n: Self) -> Self {
        self.powf(n)
    }
}

#[test]
fn test_float_cbrt() {
    assert_eq!(-8f64.cbrt(), -2f64);
    assert_eq!(8f64.cbrt(), 2f64);
    assert_eq!(0f32.cbrt(), 0f32);
}
