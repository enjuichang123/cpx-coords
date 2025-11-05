// Copyright 2025 En-Jui Chang
//
// Licensed under the Apache License, Version 2.0 <http://www.apache.org/licenses/LICENSE-2.0>
// or the MIT license <http://opensource.org/licenses/MIT>, at your option.
// This file may not be copied, modified, or distributed except according to those terms.

use core::hash::Hash;
use core::ops::{Add, Div, Mul, Neg, Sub};
use num_traits::{Float, NumCast};

#[doc(hidden)]
/// Private module to seal the `f32` or `f64`.
mod float_sealed {
    pub trait FloatSealed {}
    impl FloatSealed for f32 {}
    impl FloatSealed for f64 {}
}

/// Extension trait for providing common mathematical constants for floats.
pub trait FloatExt:
    float_sealed::FloatSealed + Float + NumCast + Copy + PartialOrd + 'static
{
    const ZERO: Self;
    const ONE: Self;
    const NEG_ONE: Self;
    const TWO: Self;
    const SQRT_2: Self;
    const FRAC_1_SQRT_2: Self;
    const PI: Self;
    const FRAC_1_PI: Self;
    const FRAC_PI_2: Self;
    const FRAC_2_PI: Self;
    const TAU: Self;
    const E: Self;
    const FRAC_1_E: Self;
    const FRAC_2_SQRT_PI: Self;
    const NEG_FRAC_PI_2: Self;
    const FRAC_PI_3: Self;
    const FRAC_PI_4: Self;
    const FRAC_PI_6: Self;
    const FRAC_PI_8: Self;
    const LN_2: Self;
    const LN_10: Self;
    const LOG2_10: Self;
    const LOG2_E: Self;
    const LOG10_2: Self;
    const LOG10_E: Self;
    const EXP_PI: Self;
    const EXP_NEG_PI: Self;
    const EXP_FRAC_PI_2: Self;
    const EXP_NEG_FRAC_PI_2: Self;
    const EPSILON: Self;
    /// Returns the threshold below which two floats are considered equal.
    ///
    /// # Examples
    ///
    /// ```
    /// use cpx_coords::FloatExt;
    ///
    /// assert_eq!(f32::THRESHOLD,10.0 * f32::EPSILON);
    /// assert_eq!(f64::THRESHOLD,10.0 * f64::EPSILON);
    /// ```
    const THRESHOLD: Self;
    /// Returns the `ln(threshold)`.
    ///
    /// # Examples
    ///
    /// ```
    /// use cpx_coords::FloatExt;
    ///
    /// assert_eq!(f32::LN_THRESHOLD,-13.6398001_f32);
    /// assert_eq!(f64::LN_THRESHOLD,-67.4821365922_f64);
    /// ```
    const LN_THRESHOLD: Self;

    /// Approximately checks whether two `T` numbers are equal within a given threshold.
    ///
    /// # Arguments
    /// * `a` - First floating-point number.
    /// * `b` - Second floating-point number.
    ///
    /// # Returns
    /// * `true` if the absolute difference between `a` and `b` is less than or equal to `threshold`, otherwise `false`.
    ///
    /// # Examples
    /// ```
    /// use cpx_coords::FloatExt;
    ///
    /// type F1 = f32;
    /// let a1: F1 = 1.0;
    /// let b1: F1 = a1 + 0.9 * F1::THRESHOLD;
    /// assert!(a1.approx_eq(b1));
    /// let c1: F1 = a1 + 1.1 * F1::THRESHOLD;
    /// assert!(!a1.approx_eq(c1));
    ///
    /// type F2 = f64;
    /// let a2: F2 = 1.0;
    /// let b2: F2 = a2 - 1.0 * F2::THRESHOLD;
    /// assert!(a2.approx_eq(b2));
    /// let c2: F2 = a2 + 1.1 * F2::THRESHOLD;
    /// assert!(!a2.approx_eq(c2));
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f32;
    /// let a: F = 1.0;
    /// let b: F = 2.0;
    ///
    /// assert!((..=4).contains(&measure_cycles( || { let _ = a.approx_eq(b); }, || { let _ = a; }, batch, reps)));
    /// ```
    fn approx_eq(self, other: Self) -> bool;

    /// Wraps a phase angle into the `(-π, π]` interval.
    ///
    /// # Examples
    /// ```
    /// use cpx_coords::FloatExt;
    ///
    /// type F = f32;
    ///
    /// assert_eq!(0.0.principal_val(), 0.0);
    /// assert_eq!(F::PI.principal_val(), F::PI);
    /// assert_eq!((-F::PI).principal_val(), F::PI); // -π maps to π
    ///
    ///
    /// let x = (F::PI + 0.1).principal_val();
    /// let expected = -F::PI + 0.1;
    /// assert!((x - expected).abs() <= F::THRESHOLD);
    ///
    /// let y = (-F::PI - 0.1).principal_val();
    /// let expected = F::PI - 0.1;
    /// assert!((y - expected).abs() <= F::THRESHOLD);
    ///
    /// assert!((3.0 * F::TAU).principal_val().abs() <= F::THRESHOLD);
    /// ```
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// let val: F = 3.0;
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = val.principal_val(); }, || { let _ = val; }, batch, reps)));
    ///
    /// let near_val: F = -F::PI;
    /// // assert!((..=11).contains(&measure_cycles( || { let _ = near_val.principal_val(); }, || { let _ = near_val; }, batch, reps)));
    ///
    /// let mid_val: F = 15.0;
    /// // assert!((..=30).contains(&measure_cycles( || { let _ = mid_val.principal_val(); }, || { let _ = mid_val; }, batch, reps)));
    ///
    /// let far_val: F = F::MAX;
    /// // assert!((..=1850).contains(&measure_cycles( || { let _ = far_val.principal_val(); }, || { let _ = far_val; }, batch, reps)));
    /// ```
    #[inline(always)]
    fn principal_val(self) -> Self {
        if self > -Self::PI && self <= Self::PI {
            return self;
        } else if self == -Self::PI {
            return Self::PI;
        }
        let mut r = self % Self::TAU;
        if r <= -Self::PI {
            r = r + Self::TAU
        } else if r > Self::PI {
            r = r - Self::TAU
        }
        r
    }
    /// Returns the maximum of this floating type.
    const MAX: Self;
    /// Returns the minimum of this floating type.
    const MIN: Self;
    /// Returns the `ln(T::MAX)` of this floating type.
    const LN_MAX: Self;
    const INFINITY: Self;
    const NEG_INFINITY: Self;
}

impl FloatExt for f32 {
    const ZERO: Self = 0.0_f32;
    const ONE: Self = 1.0_f32;
    const NEG_ONE: Self = -1.0_f32;
    const TWO: Self = 2.0_f32;
    const SQRT_2: Self = core::f32::consts::SQRT_2;
    const FRAC_1_SQRT_2: Self = core::f32::consts::FRAC_1_SQRT_2;
    const PI: Self = core::f32::consts::PI;
    const FRAC_1_PI: Self = core::f32::consts::FRAC_1_PI;
    const FRAC_PI_2: Self = core::f32::consts::FRAC_PI_2;
    const FRAC_2_PI: Self = core::f32::consts::FRAC_2_PI;
    const TAU: Self = core::f32::consts::TAU;
    const E: Self = core::f32::consts::E;
    const FRAC_1_E: Self = core::f32::consts::E.recip();
    const FRAC_2_SQRT_PI: Self = core::f32::consts::FRAC_2_SQRT_PI;
    const NEG_FRAC_PI_2: Self = -core::f32::consts::FRAC_PI_2;
    const FRAC_PI_3: Self = core::f32::consts::FRAC_PI_3;
    const FRAC_PI_4: Self = core::f32::consts::FRAC_PI_4;
    const FRAC_PI_6: Self = core::f32::consts::FRAC_PI_6;
    const FRAC_PI_8: Self = core::f32::consts::FRAC_PI_8;
    const LN_2: Self = core::f32::consts::LN_2;
    const LN_10: Self = core::f32::consts::LN_10;
    const LOG2_10: Self = core::f32::consts::LOG2_10;
    const LOG2_E: Self = core::f32::consts::LOG2_E;
    const LOG10_2: Self = core::f32::consts::LOG10_2;
    const LOG10_E: Self = core::f32::consts::LOG10_E;
    const EXP_PI: Self = 23.140696_f32;
    const EXP_NEG_PI: Self = 0.043213915_f32;
    const EXP_FRAC_PI_2: Self = 4.8104777_f32;
    const EXP_NEG_FRAC_PI_2: Self = 0.20787957_f32;
    const EPSILON: Self = f32::EPSILON;
    const THRESHOLD: Self = f32::EPSILON * 10.0;
    const LN_THRESHOLD: Self = -13.639_8_f32;
    fn approx_eq(self, other: f32) -> bool {
        let val = self - other;
        if val > Self::THRESHOLD {
            false
        } else {
            val >= -Self::THRESHOLD
        }
        //(self - other).abs() <= Self::THRESHOLD
    }
    const MAX: Self = f32::MAX;
    const MIN: Self = f32::MIN;
    const LN_MAX: Self = 88.72284_f32;
    const INFINITY: Self = f32::INFINITY;
    const NEG_INFINITY: Self = f32::NEG_INFINITY;
}

impl FloatExt for f64 {
    const ZERO: Self = 0.0_f64;
    const ONE: Self = 1.0_f64;
    const NEG_ONE: Self = -1.0_f64;
    const TWO: Self = 2.0_f64;
    const SQRT_2: Self = core::f64::consts::SQRT_2;
    const FRAC_1_SQRT_2: Self = core::f64::consts::FRAC_1_SQRT_2;
    const PI: Self = core::f64::consts::PI;
    const FRAC_1_PI: Self = core::f64::consts::FRAC_1_PI;
    const FRAC_PI_2: Self = core::f64::consts::FRAC_PI_2;
    const FRAC_2_PI: Self = core::f64::consts::FRAC_2_PI;
    const TAU: Self = core::f64::consts::TAU;
    const E: Self = core::f64::consts::E;
    const FRAC_1_E: Self = core::f64::consts::E.recip();
    const FRAC_2_SQRT_PI: Self = core::f64::consts::FRAC_2_SQRT_PI;
    const NEG_FRAC_PI_2: Self = -core::f64::consts::FRAC_PI_2;
    const FRAC_PI_3: Self = core::f64::consts::FRAC_PI_3;
    const FRAC_PI_4: Self = core::f64::consts::FRAC_PI_4;
    const FRAC_PI_6: Self = core::f64::consts::FRAC_PI_6;
    const FRAC_PI_8: Self = core::f64::consts::FRAC_PI_8;
    const LN_2: Self = core::f64::consts::LN_2;
    const LN_10: Self = core::f64::consts::LN_10;
    const LOG2_10: Self = core::f64::consts::LOG2_10;
    const LOG2_E: Self = core::f64::consts::LOG2_E;
    const LOG10_2: Self = core::f64::consts::LOG10_2;
    const LOG10_E: Self = core::f64::consts::LOG10_E;
    const EXP_PI: Self = 23.140692632779267_f64;
    const EXP_NEG_PI: Self = 0.04321391826377226_f64;
    const EXP_FRAC_PI_2: Self = 4.810477380965351_f64;
    const EXP_NEG_FRAC_PI_2: Self = 0.20787957635076193_f64;
    const EPSILON: Self = f64::EPSILON;
    const THRESHOLD: Self = f64::EPSILON * 10.0;
    const LN_THRESHOLD: Self = -67.4821365922_f64;
    fn approx_eq(self, other: f64) -> bool {
        let val = self - other;
        if val > Self::THRESHOLD {
            false
        } else {
            val >= -Self::THRESHOLD
        }
        //(self - other).abs() <= Self::THRESHOLD
    }
    const MAX: Self = f64::MAX;
    const MIN: Self = f64::MIN;
    const LN_MAX: Self = 709.782712893384_f64;
    const INFINITY: Self = f64::INFINITY;
    const NEG_INFINITY: Self = f64::NEG_INFINITY;
}

/// Enum representing complex numbers in various coordinate systems.
#[derive(Debug, Clone, Copy, Default, Eq)]
pub enum Cpx<T>
where
    T: FloatExt,
{
    #[default]
    /// Represents the additive identity, `0.0`.
    Zero,
    /// Represents the multiplicative identity, `1.0`.
    One,
    /// Represents `-1.0`.
    NegOne,
    /// Represents the imaginary unit, `j * 1.0`, which is one of the square root of `-1.0`.
    J,
    /// Represents `j * (-1.0)`.
    NegJ,
    /// Represents a real floating number `re`.
    Re { re: T },
    /// Represents a imaginary floating number `j * im`.
    Im { im: T },
    /// Represents a complex number with the unit radius and a real floating phase angle `ph` in `(-π, π]`, i.e. `cos(ph) + j * sin(ph)`.
    Ph { ph: T },
    /// Represents a complex number in Cartesian coordinates `re + j * im`.
    ReIm { re: T, im: T },
    /// Represents a complex number in polar coordinates: `rad * exp(j * ph)`.
    Pl { rad: T, ph: T },
    /// Represents a complex number `z = exp(re + j * im)` in logarithmic form: `ln(z) = re + j * im`.
    Ln { re: T, im: T },
}

/// Canonicalize an axis expression into symbolic constants if possible.
#[macro_export]
macro_rules! canonicalize_axis {
    ($variant:ident, $val:ident, $one:ident, $neg_one:ident) => {
        match $val {
            _ if $val == T::ONE => $one,
            _ if $val == T::ZERO => Cpx::Zero,
            _ if $val == T::NEG_ONE => $neg_one,
            _ => $variant { $val },
        }
    };
}

/// Canonicalizes an axis-aligned value into a symbolic constant when possible.
///
/// This macro converts `$val` into one of the predefined constants (`$one` or `$neg_one`)
/// if it matches `T::ONE` or `T::NEG_ONE`; otherwise, it constructs the variant `$variant { $val }`.
///
/// Assumes that `$val == T::ZERO` has already been checked and excluded by the caller.
#[macro_export]
macro_rules! canonicalize_axis_skip_zero {
    ($variant:ident, $val:ident, $one:ident, $neg_one:ident) => {
        match $val {
            _ if $val == T::ONE => $one,
            _ if $val == T::NEG_ONE => $neg_one,
            _ => $variant { $val },
        }
    };
}

/// Canonicalize a phase expression into symbolic constants if possible.
#[macro_export]
macro_rules! canonicalize_phase {
    ($ph_expr:expr) => {
        match ($ph_expr).principal_val() {
            ph if ph == T::PI => Cpx::NegOne,
            ph if ph == T::FRAC_PI_2 => Cpx::J,
            ph if ph == T::ZERO => Cpx::One,
            ph if ph == T::NEG_FRAC_PI_2 => Cpx::NegJ,
            ph => Cpx::Ph { ph },
        }
    };
}
#[macro_export]
macro_rules! canonicalize_polar {
    ($rad:expr, $ph:expr) => {
        match $rad {
            r if r == T::ZERO => Cpx::Zero,
            r if r == T::ONE => canonicalize_phase!($ph),
            r if r > T::ZERO => match $ph.principal_val() {
                ph if ph == T::PI => Cpx::Re { re: -r },
                ph if ph == T::FRAC_PI_2 => Cpx::Im { im: r },
                ph if ph == T::ZERO => Cpx::Re { re: r },
                ph if ph == T::NEG_FRAC_PI_2 => Cpx::Im { im: -r },
                ph => Cpx::Pl { rad: r, ph },
            },
            r => match ($ph + T::PI).principal_val() {
                ph if ph == T::PI => Cpx::Re { re: r },
                ph if ph == T::FRAC_PI_2 => Cpx::Im { im: -r },
                ph if ph == T::ZERO => Cpx::Re { re: -r },
                ph if ph == T::NEG_FRAC_PI_2 => Cpx::Im { im: r },
                ph => Cpx::Pl { rad: -r, ph },
            },
        }
    };
}

impl<T: FloatExt> Cpx<T> {
    /// Returns the real part of the complex number as T.
    ///
    /// # Example
    /// ```rust
    /// use cpx_coords::{Cpx,FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f32;
    /// let (val,val_2): (F,F) = (3.0,4.0);
    ///
    /// assert_eq!((Zero::<F>.re(), One::<F>.re(), NegOne::<F>.re(), J::<F>.re(), NegJ::<F>.re(), Im { im: val }.re()), (0.0, 1.0, -1.0, 0.0, 0.0, 0.0));
    /// assert_eq!((Re { re: val }.re(), ReIm { re: val, im: val_2 }.re()),(val,val));
    /// assert_eq!((Ph { ph: val_2 }.re(), Pl { rad: val, ph: val_2 }.re(), Ln { re: val, im: val_2 }.re()),(val_2.cos(), val * val_2.cos(), val.exp() * val_2.cos()));
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (1.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_3;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = zero.re(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = re.re(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = re_im.re(); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=42).contains(&measure_cycles( || { let _ = ph.re(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=44).contains(&measure_cycles( || { let _ = pl.re(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=66).contains(&measure_cycles( || { let _ = ln.re(); }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn re(&self) -> T {
        use Cpx::*;
        match self {
            Zero | J | NegJ | Im { .. } => T::ZERO,
            One => T::ONE,
            NegOne => T::NEG_ONE,
            Re { re } | ReIm { re, .. } => *re,
            Ph { ph } => match *ph {
                x if (x == T::FRAC_PI_2) || (x == T::NEG_FRAC_PI_2) => T::ZERO, // (x == T::FRAC_PI_2) || (x == T::NEG_FRAC_PI_2) => T::ZERO,
                x if x == T::ZERO => T::ONE,
                x if x == T::PI => T::NEG_ONE,
                x => x.cos(),
            },
            Pl { rad, ph } => match *ph {
                x if (x == T::FRAC_PI_2) || (x == T::NEG_FRAC_PI_2) => T::ZERO,
                x if x == T::ZERO => *rad,
                x if x == T::PI => -*rad,
                x => x.cos() * *rad,
            },
            Ln { re, im } => match *re {
                r if (r < T::LN_THRESHOLD)
                    || (*im == T::FRAC_PI_2)
                    || (*im == T::NEG_FRAC_PI_2) =>
                {
                    T::ZERO
                }
                r if r >= T::LN_MAX => match im.principal_val() {
                    i if (i < T::NEG_FRAC_PI_2) || (i > T::FRAC_PI_2) => T::NEG_INFINITY,
                    i if (i == T::NEG_FRAC_PI_2) || (i == T::FRAC_PI_2) => T::ZERO,
                    _ => T::INFINITY,
                },
                r => {
                    let rad = r.exp();
                    match *im {
                        i if i == T::ZERO => rad,
                        i if i == T::PI => -rad,
                        i => i.cos() * rad,
                    }
                }
            },
        }
    }

    /// Returns the imaginary part of the complex number as T.
    ///
    /// # Example
    /// ```rust
    /// use cpx_coords::{Cpx,FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f64;
    /// let val_1: F = 3.0;
    /// let val_2: F = 4.0;
    ///
    /// assert_eq!((Zero::<F>.im(), One::<F>.im(), NegOne::<F>.im(), J::<F>.im(), NegJ::<F>.im(), Im { im: val_1 }.im()), (0.0, 0.0, 0.0, 1.0, -1.0, val_1));
    /// assert_eq!((Re { re: val_1 }.im(), ReIm { re: val_1, im: val_2 }.im()),(0.0,val_2));
    /// assert_eq!((Ph { ph: val_2 }.im(), Pl { rad: val_1, ph: val_2 }.im(), Ln { re: val_1, im: val_2 }.im()),(val_2.sin(), val_1 * val_2.sin(), val_1.exp() * val_2.sin()));
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_3;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = zero.im(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = im.im(); }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = re_im.im(); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=38).contains(&measure_cycles( || { let _ = ph.im(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=42).contains(&measure_cycles( || { let _ = pl.im(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=65).contains(&measure_cycles( || { let _ = ln.im(); }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn im(&self) -> T {
        use Cpx::*;
        match self {
            Zero | One | NegOne | Re { .. } => T::ZERO,
            J => T::ONE,
            NegJ => T::NEG_ONE,
            Im { im } | ReIm { im, .. } => *im,
            Ph { ph } => match *ph {
                x if (x == T::ZERO) || (x == T::PI) => T::ZERO,
                x if x == T::FRAC_PI_2 => T::ONE,
                x if x == T::NEG_FRAC_PI_2 => T::NEG_ONE,
                x => x.sin(),
            },
            Pl { rad, ph } => match *ph {
                x if (x == T::ZERO) || (x == T::PI) => T::ZERO,
                x if x == T::FRAC_PI_2 => *rad,
                x if x == T::NEG_FRAC_PI_2 => -*rad,
                x => x.sin() * *rad,
            },
            Ln { re, im } => match *re {
                r if (r < T::LN_THRESHOLD) || (*im == T::ZERO) || (*im == T::PI) => T::ZERO,
                r if r >= T::LN_MAX => match im.principal_val() {
                    i if (i == T::ZERO) || (i == T::PI) => T::ZERO,
                    i if i < T::ZERO => T::NEG_INFINITY,
                    _ => T::INFINITY,
                },
                r => {
                    let rad = r.exp();
                    match *im {
                        i if i == T::FRAC_PI_2 => rad,
                        i if i == T::NEG_FRAC_PI_2 => -rad,
                        i => i.sin() * rad,
                    }
                }
            },
        }
    }

    /// Returns the real and imaginary components of the complex number as a tuple `(re, im)`.
    ///
    /// Both components are of type `T`. This method is equivalent to calling `.re()` and `.im()` separately,
    /// but returns the result in a single call for performance benefits.
    ///
    ///
    /// # Example
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f64;
    /// let val_1: F = 3.0;
    /// let val_2: F = 4.0;
    ///
    /// assert_eq!((Zero::<F>.re_im(), One::<F>.re_im(), NegOne::<F>.re_im(), J::<F>.re_im(), NegJ::<F>.re_im()),((0.0, 0.0), (1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)));
    /// assert_eq!((Re { re: val_1 }.re_im(), Im { im: val_2 }.re_im(), ReIm { re: val_1, im: val_2 }.re_im()),((val_1, 0.0), (0.0, val_2), (val_1, val_2)));
    ///
    /// let ph = F::FRAC_PI_3;
    /// assert_eq!((Ph { ph }.re_im(), Ln { re: val_1, im: ph }.re_im(), Pl { rad: val_1, ph }.re_im()),((ph.cos(),ph.sin()), (val_1.exp() * ph.cos(), val_1.exp() * ph.sin()), (val_1 * ph.cos(), val_1 * ph.sin())));
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // Improvement compare to use individual `.re()` and `.im()`.
    ///
    /// // Minor improvement:
    /// // assert!((1..).contains(&measure_cycles( || { let _ = (zero.re(), zero.im()); }, || { let _ = zero.re_im(); }, batch, reps)));
    /// // assert!((2..).contains(&measure_cycles( || { let _ = (re.re(), re.im()); }, || { let _ = re.re_im(); }, batch, reps)));
    /// // assert!((3..).contains(&measure_cycles( || { let _ = (re_im.re(), re_im.im()); }, || { let _ = re_im.re_im(); }, batch, reps)));
    ///
    /// // Moderate improvement:
    /// // Computing `.cos()` and `.sin()` within the same match arm avoids redundant branching and separate match evaluations.
    /// // assert!((18..).contains(&measure_cycles( || { let _ = (ph.re(), ph.im()); }, || { let _ = ph.re_im(); }, batch, reps)));
    /// // assert!((18..).contains(&measure_cycles( || { let _ = (pl.re(), pl.im()); }, || { let _ = pl.re_im(); }, batch, reps)));
    ///
    /// // Higher improvement:
    /// // Same reason as above, with the additional benefit of avoiding two `.exp()` computations.
    /// // assert!((27..).contains(&measure_cycles( || { let _ = (ln.re(), ln.im()); }, || { let _ = ln.re_im(); }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn re_im(self) -> (T, T) {
        use Cpx::*;
        match self {
            Zero => (T::ZERO, T::ZERO),
            One => (T::ONE, T::ZERO),
            NegOne => (T::NEG_ONE, T::ZERO),
            J => (T::ZERO, T::ONE),
            NegJ => (T::ZERO, T::NEG_ONE),
            Re { re } => (re, T::ZERO),
            Im { im } => (T::ZERO, im),
            ReIm { re, im } => (re, im),
            //Ph { ph } | Pl { ph, .. } | Ln { im: ph, .. } => {
            //    let (cos, sin) = (ph.cos(), ph.sin());
            //    match self {
            //        Ph { .. } => (cos, sin),
            //        Pl { rad, .. } => (cos * rad, sin * rad),
            //        Ln { re, .. } => {
            //            let rad = re.exp();
            //            (cos * rad, sin * rad)
            //        }
            //        _ => unreachable!(),
            //    }
            //}
            // This expanded version outperforms the simpler implementation. For special cases, early returns
            // eliminate unnecessary function calls, reducing overhead significantly.
            // Even in the general case, where `.cos()`, `.sin()`, and `.exp()` are computed, the performance
            // improves (e.g., ~75 CPU cycles down to ~45).
            // The main difference is avoiding the nested match, which itself costs only 1–2 cycles, so the
            // speedup is surprising. Likely, the compiler achieves better instruction-level parallelism (ILP)
            // and branch prediction with this structure.
            // Further investigation using LLVM IR or assembly analysis might confirm this hypothesis.
            Ph { ph } => match ph {
                x if x == T::ZERO => (T::ONE, T::ZERO),
                x if x == T::PI => (T::NEG_ONE, T::ZERO),
                x if x == T::FRAC_PI_2 => (T::ZERO, T::ONE),
                x if x == T::NEG_FRAC_PI_2 => (T::ZERO, T::NEG_ONE),
                _ => (ph.cos(), ph.sin()),
            },
            Pl { rad, ph } => match ph {
                x if x == T::ZERO => (rad, T::ZERO),
                x if x == T::PI => (-rad, T::ZERO),
                x if x == T::FRAC_PI_2 => (T::ZERO, rad),
                x if x == T::NEG_FRAC_PI_2 => (T::ZERO, -rad),
                _ => (rad * ph.cos(), rad * ph.sin()),
            },
            Ln { re, im } => match re {
                r if (r < T::LN_THRESHOLD) || (im == T::ZERO) || (im == T::PI) => {
                    (T::ZERO, T::ZERO)
                }
                r if r >= T::LN_MAX => match im.principal_val() {
                    i if i == T::ZERO => (T::INFINITY, T::ZERO),
                    i if i == T::PI => (T::NEG_INFINITY, T::ZERO),
                    i if i == T::FRAC_PI_2 => (T::ZERO, T::INFINITY),
                    i if i == T::NEG_FRAC_PI_2 => (T::ZERO, T::NEG_INFINITY),
                    i if (i > -T::PI) && (i < T::NEG_FRAC_PI_2) => {
                        (T::NEG_INFINITY, T::NEG_INFINITY)
                    }
                    i if (i > T::NEG_FRAC_PI_2) && (i < T::ZERO) => (T::INFINITY, T::NEG_INFINITY),
                    i if (i > T::ZERO) && (i < T::FRAC_PI_2) => (T::INFINITY, T::INFINITY),
                    _ => (T::NEG_INFINITY, T::INFINITY),
                },
                r => {
                    let rad = r.exp();
                    match im {
                        i if i == T::ZERO => (rad, T::ZERO),
                        i if i == T::PI => (-rad, T::ZERO),
                        i if i == T::FRAC_PI_2 => (T::ZERO, rad),
                        i if i == T::NEG_FRAC_PI_2 => (T::ZERO, -rad),
                        i => (i.cos() * rad, i.sin() * rad),
                    }
                }
            },
        }
    }

    /// Converts this complex number into its canonical rectangular (real–imaginary) form.
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = zero.to_re_im(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = re.to_re_im(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = re_im.to_re_im(); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=58).contains(&measure_cycles( || { let _ = ph.to_re_im(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=63).contains(&measure_cycles( || { let _ = pl.to_re_im(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=95).contains(&measure_cycles( || { let _ = ln.to_re_im(); }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn to_re_im(self) -> Self {
        use Cpx::*;
        let (re, im) = self.re_im();
        ReIm { re, im }
    }

    /// Canonicalizes and simplifies a `Cpx<T>` to its rectangular (real–imaginary) form if possible.
    /// This method normalizes the representation of the complex number to one of the canonical forms:
    /// - `Zero`, `One`, `NegOne`, `J`, `NegJ`
    /// - `Re { re }`
    /// - `Im { im }`
    /// - `ReIm { re, im }`
    ///
    /// # Performance Considerations
    /// - [`Cpx::canonicalize_re_im()`] is generally **more expensive** (~10 CPU cycles) than [`Cpx::to_re_im()`]
    ///   because it performs extra checks and canonicalization.
    /// - For simple equality checks, [`Cpx::to_re_im()`] is usually sufficient.
    /// - For **addition**, [`Cpx::to_re_im()`] typically provides normalized operands and is preferred in most cases.
    /// - Use [`Cpx::canonicalize_re_im()`] only when canonical variants (`Zero`, `One`, `NegOne`, `J`, `NegJ`,
    ///   `Re { re }`, `Im { im }`) are highly desirable.
    ///   - Example: when summing many complex numbers, grouping pure real and imaginary parts first
    ///     allows extremely fast aggregation, followed by combining the remaining components.
    ///   - Another example: in non-performance-critical cases, such as rendering complex numbers as
    ///     LaTeX strings, canonicalization can simplify formatting.
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = zero.canonicalize_re_im(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=14).contains(&measure_cycles( || { let _ = re.canonicalize_re_im(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = re_im.canonicalize_re_im(); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=71).contains(&measure_cycles( || { let _ = ph.canonicalize_re_im(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=79).contains(&measure_cycles( || { let _ = pl.canonicalize_re_im(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=104).contains(&measure_cycles( || { let _ = ln.canonicalize_re_im(); }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn canonicalize_re_im(self) -> Self {
        use Cpx::*;
        match self {
            Zero | One | NegOne | J | NegJ => self,
            Re { re } => canonicalize_axis!(Re, re, One, NegOne),
            Im { im } => canonicalize_axis!(Im, im, J, NegJ),
            ReIm { re, im } => {
                if im == T::ZERO {
                    canonicalize_axis!(Re, re, One, NegOne)
                } else if re == T::ZERO {
                    canonicalize_axis_skip_zero!(Im, im, J, NegJ)
                } else {
                    self
                }
            }
            _ => {
                let (re, im) = self.re_im();
                if im == T::ZERO {
                    canonicalize_axis!(Re, re, One, NegOne)
                } else if re == T::ZERO {
                    canonicalize_axis_skip_zero!(Im, im, J, NegJ)
                } else {
                    self
                }
            }
        }
    }

    /// Returns the radius (magnitude) of the complex number as T.
    ///
    /// For rectangular coordinates (`ReIm { re, im }`), this computes
    /// `sqrt(re² + im²)` using [`f32::hypot`] or [`f64::hypot`], which minimizes intermediate
    /// overflow/underflow. However, if both `re` and `im` are near `T::MAX`,
    /// the result will exceed `T::MAX` and return positive infinity (`∞`).
    ///
    /// # Example
    /// ```rust
    /// use cpx_coords::{Cpx,FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f64;
    ///
    /// let ph = F::FRAC_PI_3;
    /// assert_eq!((Zero::<F>.rad(), One::<F>.rad(), NegOne::<F>.rad(), J::<F>.rad(), NegJ::<F>.rad(), Ph::<F> { ph }.rad()),(0.0, 1.0, 1.0, 1.0, 1.0, 1.0));
    ///
    /// let (val, val_2): (F,F) = (-3.0, 4.0);
    /// assert_eq!((Re { re: val }.rad(), Im { im: val }.rad(), ReIm { re: val, im: val_2 }.rad()),(val.abs(), val.abs(), val.hypot(val_2)));
    ///
    /// let rad: F = 2.0;
    /// assert_eq!((Ln { re: rad.ln(), im: ph }.rad(), Pl { rad, ph }.rad()),(rad, rad));
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (mut zero, mut one, mut neg_one, mut j, mut neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (mut re, mut im, mut re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (mut ph, mut pl, mut ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = zero.rad(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = ph.rad(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = pl.rad(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = re.rad(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=17).contains(&measure_cycles( || { let _ = re_im.rad(); }, || { let _ = re_im; }, batch, reps)));
    ///  assert!((..=23).contains(&measure_cycles( || { let _ = ln.rad(); }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn rad(&self) -> T {
        use Cpx::*;
        match self {
            Zero => T::ZERO,
            One | NegOne | J | NegJ | Ph { .. } => T::ONE,
            Pl { rad, .. } => *rad,
            Re { re } | Im { im: re } => {
                if *re >= T::ZERO {
                    *re
                } else {
                    -*re
                }
            }
            ReIm { re, im } => re.hypot(*im),
            Ln { re, .. } => match *re {
                x if x < T::LN_THRESHOLD => T::ZERO,
                x if x > T::LN_MAX => T::INFINITY,
                x => x.exp(),
            },
        }
    }

    /// Returns the phase angle of the complex number as `T`.
    ///
    /// # Example
    /// ```rust
    /// use cpx_coords::{Cpx,FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f64;
    ///
    /// //assert_eq!((Zero::<F>.ph(), One::<F>.ph(), NegOne::<F>.ph(), J::<F>.ph(), NegJ::<F>.ph()), (0.0, 0.0, F::PI, F::FRAC_PI_2, F::NEG_FRAC_PI_2));
    ///
    /// let pos_val: F = 3.0;
    /// //assert_eq!((Re { re: pos_val }.ph(), Re { re: -pos_val }.ph(), Im { im: pos_val }.ph(), Im { im: -pos_val }.ph()), (0.0, F::PI, F::FRAC_PI_2, F::NEG_FRAC_PI_2));
    ///
    /// let ph: F = F::FRAC_PI_6.principal_val();
    /// //assert_eq!((Ph { ph }.ph(), Pl { rad: pos_val, ph }.ph(), Ln { re: pos_val, im: ph }.ph()), (ph, ph, ph));
    ///
    /// let (re, im): (F, F) = (3.0, 4.0);
    /// //assert_eq!(ReIm { re, im }.ph(), im.atan2(re));
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (mut zero, mut one, mut neg_one, mut j, mut neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (mut re, mut im, mut re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_2;
    /// let (mut ph, mut pl, mut ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = zero.ph(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = re.ph(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = ph.ph(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=7).contains(&measure_cycles( || { let _ = pl.ph(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=7).contains(&measure_cycles( || { let _ = ln.ph(); }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=33).contains(&measure_cycles( || { let _ = re_im.ph(); }, || { let _ = re_im; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn ph(&self) -> T {
        use Cpx::*;
        match self {
            Zero | One => T::ZERO,
            J => T::FRAC_PI_2,
            NegOne => T::PI,
            NegJ => T::NEG_FRAC_PI_2,
            Re { re } => match *re {
                r if r >= T::ZERO => T::ZERO,
                _ => T::PI,
            },
            Im { im } => match *im {
                i if i < T::ZERO => T::NEG_FRAC_PI_2,
                i if i > T::ZERO => T::FRAC_PI_2,
                _ => T::ZERO,
            },
            Ph { ph } => ph.principal_val(),
            Pl { rad, ph } => match *rad {
                rad if rad > T::ZERO => ph.principal_val(),
                _ => T::ZERO,
            },
            Ln { re, im: ph } => match *re {
                r if r >= T::LN_THRESHOLD => ph.principal_val(),
                _ => T::ZERO,
            },
            ReIm { re, im } => im.atan2(*re),
        }
    }

    /// Returns the polar coordinates `(r, φ)` of the complex number, where:
    ///
    /// * `r` is the radius (magnitude)
    /// * `φ` is the phase angle in radians
    ///
    /// Both values are returned as type `T`.
    ///
    /// For rectangular coordinates (`ReIm { re, im }`), the radius is computed as
    /// `sqrt(re² + im²)` using [`f32::hypot`] or [`f64::hypot`], which minimizes intermediate
    /// overflow and underflow. However, if both `re` and `im` are close to
    /// [`FloatExt::MAX`], the result will exceed the representable range and return
    /// positive infinity (`∞`).
    ///
    /// The phase is normalized using the [`FloatExt::principal_val()`] form to ensure a
    /// well-defined principal value in the range `(-π, π]`.
    ///
    /// # Example
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f64;
    ///
    /// assert_eq!((Zero::<F>.rad_ph(), One::<F>.rad_ph(), NegOne::<F>.rad_ph(), J::<F>.rad_ph(), NegJ::<F>.rad_ph()), ((0.0, 0.0), (1.0, 0.0), (1.0, F::PI), (1.0, F::FRAC_PI_2), (1.0, F::NEG_FRAC_PI_2)));
    ///
    /// let pos_val: F = 3.0;
    /// assert_eq!((Re { re: pos_val }.rad_ph(), Re { re: -pos_val }.rad_ph(), Im { im: pos_val }.rad_ph(), Im { im: -pos_val }.rad_ph()), ((pos_val, 0.0), (pos_val, F::PI), (pos_val, F::FRAC_PI_2), (pos_val, F::NEG_FRAC_PI_2)));
    ///
    /// let (rad, ph): (F, F) = (4.0, F::FRAC_PI_6.principal_val());
    /// assert_eq!((Ph { ph }.rad_ph(), Pl { rad, ph }.rad_ph(), Ln { re: rad.ln(), im: ph }.rad_ph()), ((1.0, ph), (rad, ph), (rad, ph)));
    ///
    /// let (re, im): (F, F) = (3.0, 4.0);
    /// assert_eq!(ReIm { re, im }.rad_ph(), (re.hypot(im), im.atan2(re)));
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // The CPU cycles improved comparing to separate calls.
    /// // assert!((2..).contains(&measure_cycles( || { let _ = (zero.rad(), zero.ph()); }, || { let _ = zero.rad_ph(); }, batch, reps)));
    /// // assert!((1..).contains(&measure_cycles( || { let _ = (ph.rad(), ph.ph()); }, || { let _ = ph.rad_ph(); }, batch, reps)));
    /// // assert!((0..).contains(&measure_cycles( || { let _ = (pl.rad(), pl.ph()); }, || { let _ = pl.rad_ph(); }, batch, reps)));
    /// // assert!((2..).contains(&measure_cycles( || { let _ = (ln.rad(), ln.ph()); }, || { let _ = ln.rad_ph(); }, batch, reps)));
    /// // assert!((3..).contains(&measure_cycles( || { let _ = (re.rad(), re.ph()); }, || { let _ = re.rad_ph(); }, batch, reps)));
    /// // assert!((3..).contains(&measure_cycles( || { let _ = (re_im.rad(), re_im.ph()); }, || { let _ = re_im.rad_ph(); }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn rad_ph(self) -> (T, T) {
        use Cpx::*;
        match self {
            Zero => (T::ZERO, T::ZERO),
            One => (T::ONE, T::ZERO),
            J => (T::ONE, T::FRAC_PI_2),
            NegOne => (T::ONE, T::PI),
            NegJ => (T::ONE, T::NEG_FRAC_PI_2),
            Re { re } => match re {
                r if r < T::ZERO => (-r, T::PI),
                r => (r, T::ZERO),
            },
            Im { im } => match im {
                i if i < T::ZERO => (-i, T::NEG_FRAC_PI_2),
                i if i > T::ZERO => (i, T::FRAC_PI_2),
                _ => (T::ZERO, T::ZERO),
            },
            Ph { ph } => (T::ONE, ph.principal_val()),
            Pl { rad, ph } => match rad {
                _ if rad > T::ZERO => (rad, ph.principal_val()),
                _ => (T::ZERO, T::ZERO),
            },
            Ln { re, im } => match re {
                r if r < T::LN_THRESHOLD => (T::ZERO, T::ZERO),
                r if r > T::LN_MAX => (T::INFINITY, im.principal_val()),
                r => (r.exp(), im.principal_val()),
            },
            ReIm { re, im } => (re.hypot(im), im.atan2(re)),
        }
    }

    /// Converts this complex number into its canonical polar (radius–phase) form.
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = zero.to_polar(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = re.to_polar(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=9).contains(&measure_cycles( || { let _ = ph.to_polar(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=9).contains(&measure_cycles( || { let _ = pl.to_polar(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=32).contains(&measure_cycles( || { let _ = ln.to_polar(); }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=54).contains(&measure_cycles( || { let _ = re_im.to_polar(); }, || { let _ = re_im; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn to_polar(self) -> Self {
        use Cpx::*;
        // Using `rad_ph()` in a single call is about 3 CPU cycles faster than computing
        // `rad()` and `ph()` separately.
        // For simple variants, this yields a significant improvement (≈25–50%),
        // while for more complex variants the gain is smaller but still noticeable (≈5–10%).
        //let (rad, ph) = (self.rad(), self.ph());
        let (rad, ph) = self.rad_ph();
        Pl { rad, ph }
    }

    /// Converts this complex number into its canonical polar (radius–phase) form.
    /// This method normalizes all coordinate variants to one of the canonical forms:
    /// - `Zero`, `One`, `NegOne`, `J`, `NegJ`
    /// - `Re { re }`
    /// - `Im { im }`
    /// - `Ph { ph }`
    /// - `Pl { re, im }`
    ///
    /// # Performance Considerations
    /// - [`Cpx::canonicalize_polar()`] is generally **more expensive** (~14-54 CPU cycles) than [`Cpx::to_polar()`] because it performs additional checks
    ///   and conversions for canonicalization.
    /// - For **multiplication**, [`Cpx::to_polar()`] typically provides normalized operands and is preferred in most cases.
    /// - Use [`Cpx::canonicalize_polar()`] only when canonical variants (`Zero`, `One`, `NegOne`, `J`, `NegJ`,
    ///   `Re { re }`, `Im { im }`) are highly desirable.
    ///   - Example: when multiplying many complex numbers, grouping simple variants first
    ///     allows extremely fast aggregation, followed by combining the remaining components.
    ///   - Another example: in non-performance-critical cases, such as rendering complex numbers as
    ///     LaTeX strings, canonicalization can simplify formatting.
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = zero.canonicalize_polar(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=17).contains(&measure_cycles( || { let _ = re.canonicalize_polar(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=27).contains(&measure_cycles( || { let _ = ph.canonicalize_polar(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = pl.canonicalize_polar(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=73).contains(&measure_cycles( || { let _ = ln.canonicalize_polar(); }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=108).contains(&measure_cycles( || { let _ = re_im.canonicalize_polar(); }, || { let _ = re_im; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn canonicalize_polar(self) -> Self {
        use Cpx::*;
        match self {
            Zero | One | NegOne | J | NegJ => self,
            Re { re } => canonicalize_axis!(Re, re, One, NegOne),
            Im { im } => canonicalize_axis!(Im, im, J, NegJ),
            Ph { ph } => canonicalize_phase!(ph),
            Pl { rad, ph } => canonicalize_polar!(rad, ph),
            _ => {
                let (rad, ph) = self.rad_ph();
                canonicalize_polar!(rad, ph)
            }
        }
    }

    /// Canonicalizes and simplifies a `Cpx<T>` form if possible.
    /// This method normalizes the representation of the complex number to one of the canonical forms:
    /// - `Zero`, `One`, `NegOne`, `J`, `NegJ`
    /// - `Re { re }`
    /// - `Im { im }`
    /// - `ReIm { re, im }`
    /// - `Pl { rad, ph }`
    ///
    /// Returns a canonicalized form of the complex number.
    ///
    /// ## Examples
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f64;
    ///
    /// assert_eq!(Re { re: 0.0 }.canonicalize_lazy(), Zero);
    /// assert_eq!(Re { re: 1.0 }.canonicalize_lazy(), One);
    /// assert_eq!(Re { re: -1.0 }.canonicalize_lazy(), NegOne);
    /// assert_eq!(Im { im: 1.0 }.canonicalize_lazy(), J);
    /// assert_eq!(Im { im: -1.0 }.canonicalize_lazy(), NegJ);
    ///
    /// assert_eq!(Ph { ph: 0.0 }.canonicalize_lazy(), One);
    /// assert_eq!(Ph { ph: F::FRAC_PI_2 }.canonicalize_lazy(), J);
    /// assert_eq!(Ph { ph: -F::PI }.canonicalize_lazy(), NegOne);
    ///
    /// let ph = F::FRAC_PI_2;
    /// assert_eq!(Pl { rad: 1.0, ph }.canonicalize_lazy(), J);
    /// assert_eq!(Pl { rad: 2.0, ph }.canonicalize_lazy(), Im { im: 2.0 });
    ///
    /// let z = Ln { re: 1.0_f64.ln(), im: 0.0 };
    /// assert_eq!(z.canonicalize_lazy(), One);
    ///
    /// let z = Ln { re: 2.0_f64.ln(), im: F::FRAC_PI_2 };
    /// assert_eq!(z.canonicalize_lazy(), Im { im: 2.0 });
    /// ```
    ///
    /// # Performance Considerations
    /// - [`Cpx::canonicalize_lazy()`] has more moderate worst-case performance than [`Cpx::canonicalize_re_im()`] or [`Cpx::canonicalize_polar()`],
    ///   as it adaptively chooses the most appropriate canonical form based on the current representation, avoiding costly
    ///   conversions (e.g., from `ReIm` to `Pl`).
    /// - Prefer [`Cpx::canonicalize_lazy()`] when canonical forms such as `Zero`, `One`, `NegOne`, `J`, `NegJ`,
    ///   `Re { re }`, or `Im { im }` are strongly desired, or when implementing new functions that benefit from
    ///   structured, predictable inputs with moderate worst-case performance.
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = zero.canonicalize_lazy(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=16).contains(&measure_cycles( || { let _ = re.canonicalize_lazy(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=12).contains(&measure_cycles( || { let _ = re_im.canonicalize_lazy(); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=25).contains(&measure_cycles( || { let _ = ph.canonicalize_lazy(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=35).contains(&measure_cycles( || { let _ = pl.canonicalize_lazy(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=58).contains(&measure_cycles( || { let _ = ln.canonicalize_lazy(); }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn canonicalize_lazy(self) -> Self {
        use Cpx::*;
        match self {
            Zero | One | NegOne | J | NegJ => self,
            Re { re } => canonicalize_axis!(Re, re, One, NegOne),
            Im { im } => canonicalize_axis!(Im, im, J, NegJ),
            ReIm { re, im } => {
                if im == T::ZERO {
                    canonicalize_axis!(Re, re, One, NegOne)
                } else if re == T::ZERO {
                    canonicalize_axis_skip_zero!(Im, im, J, NegJ)
                } else {
                    self
                }
            }
            Ph { ph } => canonicalize_phase!(ph),
            Pl { rad, ph } => canonicalize_polar!(rad, ph),
            Ln { im: ph, .. } => {
                let rad = self.rad();
                canonicalize_polar!(rad, ph)
            }
        }
    }

    /// Returns the complex conjugate of this number.
    ///
    /// # Examples
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt};
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// assert_eq!((C::Zero.conj(), C::One.conj(), C::NegOne.conj(), C::J.conj(), C::NegJ.conj()), (C::Zero, C::One, C::NegOne, C::NegJ, C::J));
    /// let val = 3.0;
    /// assert_eq!((C::Re { re: val }.conj(), C::Im { im: val }.conj()), (C::Re { re: val }, C::Im { im: -val }));
    /// let val_2 = 4.0;
    /// assert_eq!(C::ReIm { re: val, im: val_2 }.conj(), C::ReIm { re: val, im: -val_2 });
    ///
    /// let (rad, ph): (F, F) = (1.0, F::FRAC_PI_6);
    /// assert_eq!((C::Ph { ph }.conj(), C::Pl { rad, ph }.conj(), C::Ln { re: rad.ln(), im: ph }.conj()),(C::Ph { ph: -ph }, C::Pl { rad, ph: -ph }, C::Ln { re: rad.ln(), im: -ph }));
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = zero.conj(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = re.conj(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = neg_j.conj(); }, || { let _ = neg_j; }, batch, reps)));
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = im.conj(); }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = re_im.conj(); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=14).contains(&measure_cycles( || { let _ = ph.conj(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=12).contains(&measure_cycles( || { let _ = pl.conj(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=11).contains(&measure_cycles( || { let _ = ln.conj(); }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn conj(self) -> Self {
        use Cpx::*;
        match self {
            Zero | One | NegOne | Re { .. } => self,
            J => NegJ,
            NegJ => J,
            Im { im } => Im { im: -im },
            ReIm { re, im } => ReIm { re, im: -im },
            Ph { ph } => Ph {
                ph: (-ph).principal_val(),
            },
            Pl { rad, ph } => Pl {
                rad,
                ph: (-ph).principal_val(),
            },
            Ln { re, im } => Ln {
                re,
                im: (-im).principal_val(),
            },
        }
    }
    /// Returns the reciprocal of the complex number.
    ///
    /// This method serves as a helper for implementing division.
    ///
    /// # Examples
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt};
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// assert_eq!((C::One.recip(), C::NegOne.recip(), C::J.recip(), C::NegJ.recip()), (C::One, C::NegOne, C::NegJ, C::J));
    /// let val: F = 3.0;
    /// assert_eq!((C::Re { re: val }.recip(), C::Im { im: val }.recip()), (C::Re { re: val.recip() }, C::Im { im: -val.recip() }));
    ///
    /// let (rad, ph): (F, F) = (1.0, F::FRAC_PI_6);
    /// assert_eq!((C::Ph { ph }.recip(), C::Pl { rad, ph }.recip(), C::Ln { re: rad.ln(), im: ph }.recip()),(C::Ph { ph: -ph }, C::Pl { rad: rad.recip(), ph: -ph }, C::Ln { re: -rad.ln(), im: -ph }));
    ///
    /// let val_2: F = 4.0;
    /// let rad2: F = val.hypot(val_2).powi(2);
    /// assert_eq!(C::ReIm { re: val, im: val_2 }.recip(), C::ReIm { re: val / rad2, im: -val_2 / rad2 });
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (0.6 * F::LN_THRESHOLD,0.8);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = one.recip(); }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = neg_j.recip(); }, || { let _ = neg_j; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = zero.recip(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = ph.recip(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=8).contains(&measure_cycles( || { let _ = re.recip(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=12).contains(&measure_cycles( || { let _ = im.recip(); }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=9).contains(&measure_cycles( || { let _ = ln.recip(); }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=13).contains(&measure_cycles( || { let _ = pl.recip(); }, || { let _ = pl; }, batch, reps)));
    ///  assert!((..=25).contains(&measure_cycles( || { let _ = re_im.recip(); }, || { let _ = re_im; }, batch, reps)));
    /// ```
    #[inline(always)]
    pub fn recip(self) -> Self {
        use Cpx::*;
        match self {
            One | NegOne => self,
            J => NegJ,
            NegJ => J,
            Zero => Re { re: T::INFINITY },
            Re { re } => Re { re: re.recip() },
            Im { im } => Im { im: -im.recip() },
            Ph { ph } => Ph { ph: -ph },
            Ln { re, im } => Ln { re: -re, im: -im },
            Pl { rad, ph } => Pl {
                rad: rad.recip(),
                ph: -ph,
            },
            ReIm { re, im } => {
                //let rad2 = re.hypot(im).powi(2);
                let rad2 = re * re + im * im;
                ReIm {
                    re: re / rad2,
                    im: -im / rad2,
                }
            } //{
              //let rad = re.hypot(im);
              //if rad == T::ZERO {
              //    Re { re: T::INFINITY }
              //} else if rad == T::ONE {
              //    ReIm { re, im: -im }
              //} else {
              //    let rad2 = rad.powi(2);
              //    ReIm {
              //        re: re / rad2,
              //        im: -im / rad2,
              //    }
              //}
              //}
        }
    }

    /// Returns the exponential of the complex number.
    ///
    /// # Examples
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt};
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// // exp(0) = 1
    /// //assert_eq!(C::Zero.exp(), C::One);
    ///
    /// // exp(1) = e
    /// //assert_eq!(C::One.exp(), C::Re { re: F::E });
    ///
    /// // exp(-1) = 1/e
    /// //assert_eq!(C::NegOne.exp(), C::Re { re: F::FRAC_1_E });
    ///
    /// // exp(i) = phase(1)
    /// //assert_eq!(C::J.exp(), C::Ph { ph: 1.0 });
    ///
    /// // exp(-i) = phase(-1)
    /// //assert_eq!(C::NegJ.exp(), C::Ph { ph: -1.0 });
    ///
    /// // Real input
    /// let z = C::Re { re: 2.0 };
    /// //assert_eq!(z.exp(), C::Ln { re: 2.0, im: F::ZERO });
    ///
    /// // Imaginary input
    /// let z = C::Im { im: 2.0 };
    /// //assert_eq!(z.exp(), Cpx::Ph { ph: 2.0 });
    ///
    /// // exp(PL) → Ln(re, im)
    /// let z = C::Pl { rad: F::SQRT_2, ph: F::FRAC_PI_4 };
    /// let expected = C::Ln { re: 1.0, im: 1.0 };
    /// //assert!((z.exp() - expected).rad() < F::THRESHOLD);
    ///
    /// // exp(Ln) → Ln(re * cos, re * sin)
    /// let z = C::Ln { re: F::LN_2, im: F::FRAC_PI_4 };
    /// let expected = C::Ln { re: F::SQRT_2, im: F::SQRT_2 };
    /// //assert!((z.exp() - expected).rad() < F::THRESHOLD);
    ///
    /// // Identity of exp: 0.318132 + 1.33724i
    /// //let (v0, v1): (F,F) = (0.318132,1.3372368815372089);
    /// //let z: C = C::Pl { rad: 1.3745570107437075, ph: 1.3372357014306895 };
    /// let z: C = C::ReIm { re: 0.31813150520476413, im: 1.3372357014306895 };
    /// assert_eq!(z.exp(), z);
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = zero.exp(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = one.exp(); }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = j.exp(); }, || { let _ = j; }, batch, reps)));
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = re.exp(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = im.exp(); }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = re_im.exp(); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = ph.exp(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=47).contains(&measure_cycles( || { let _ = pl.exp(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=68).contains(&measure_cycles( || { let _ = ln.exp(); }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(never)]
    pub fn exp(self) -> Self {
        use Cpx::*;
        match self {
            Zero => One,
            One => Re { re: T::E },
            NegOne => Re { re: T::FRAC_1_E },
            J => Ph { ph: T::ONE },
            NegJ => Ph { ph: T::NEG_ONE },
            Re { re } => Ln { re, im: T::ZERO },
            Im { im: ph } => Ph { ph },
            ReIm { re, im } => Ln { re, im },
            Ph { ph } => Ln {
                re: ph.cos(),
                im: ph.sin(),
            },
            Pl { rad, ph } => Ln {
                re: rad * ph.cos(),
                im: rad * ph.sin(),
            },
            Ln { re, im } => {
                let rad = re.exp();
                Ln {
                    re: rad * im.cos(),
                    im: rad * im.sin(),
                }
            } // The above explicitly way for `Ph`, `Pl`, and `Ln` variants is about 20-30 CPU cycles faster than the below concise one.
              //x => {
              //    let (re, im) = x.re_im();
              //    Ln { re, im }
              //}
        }
    }

    /// Returns the logarithm of the complex number.
    ///
    /// # Examples
    /// ```rust
    /// use cpx_coords::{Cpx,FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f64;
    ///
    /// let zero: Cpx<F> = Zero.ln();
    /// assert!(matches!(zero, Re { re } if re.is_infinite() && re.is_sign_negative()));
    ///
    /// let one: Cpx<F> = One.ln();
    /// assert_eq!(one, Zero);
    ///
    /// let neg_one = NegOne.ln();
    /// assert_eq!(neg_one, Im { im: F::PI });
    ///
    /// let j = J.ln();
    /// assert_eq!(j, Im { im: F::FRAC_PI_2 });
    ///
    /// let neg_j = NegJ.ln();
    /// assert_eq!(neg_j, Im { im: F::NEG_FRAC_PI_2 });
    ///
    /// let re_val = Re { re: 2.0 }.ln();
    /// assert!(matches!(re_val, Re { re } if (re - F::LN_2).abs() < F::THRESHOLD));
    ///
    /// let neg_re = Re { re: -2.0 }.ln();
    /// assert!(matches!(neg_re, ReIm { re, im }
    ///     if (re - 2.0f64.ln()).abs() < F::THRESHOLD && (im - F::PI).abs() < F::THRESHOLD));
    ///
    /// let im_val = Im { im: 3.0 }.ln();
    /// assert!(matches!(im_val, ReIm { re, im }
    ///     if (re - 3.0f64.ln()).abs() < F::THRESHOLD && (im - F::FRAC_PI_2).abs() < F::THRESHOLD));
    ///
    /// let neg_im = Im { im: -4.0 }.ln();
    /// assert!(matches!(neg_im, ReIm { re, im }
    ///     if (re - 4.0f64.ln()).abs() < F::THRESHOLD && (im + F::FRAC_PI_2).abs() < F::THRESHOLD));
    ///
    /// let phase = Ph { ph: 1.2 }.ln();
    /// assert_eq!(phase, Im { im: 1.2 });
    ///
    /// let pl = Pl { rad: 2.0, ph: 1.1 }.ln();
    /// assert!(matches!(pl, ReIm { re, im }
    ///     if (re - 2.0f64.ln()).abs() < F::THRESHOLD && (im - 1.1).abs() < F::THRESHOLD));
    /// ```
    ///
    /// # Performance
    ///
    /// Cycle counts measured on `x86_64` with `measure_cycles`
    /// (`batch = 100_000`, `reps = 50`):
    ///
    /// | Variant       | Cycles (typical) |
    /// |---------------|------------------|
    /// | `Zero`        | ≤ 5              |
    /// | `One`         | ≤ 6              |
    /// | `J`           | ≤ 6              |
    /// | `Im`          | ≤ 15             |
    /// | `Re`          | ≤ 25             |
    /// | `ReIm`        | ≤ 40             |
    /// | `Ph`          | ≤ 65             |
    /// | `Pl`          | ≤ 85             |
    /// | `Ln`          | ≤ 110            |
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=4).contains(&measure_cycles( || { let _ = one.ln(); }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = j.ln(); }, || { let _ = j; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = zero.ln(); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = ph.ln(); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = ln.ln(); }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=22).contains(&measure_cycles( || { let _ = pl.ln(); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=22).contains(&measure_cycles( || { let _ = re.ln(); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=22).contains(&measure_cycles( || { let _ = im.ln(); }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=85).contains(&measure_cycles( || { let _ = re_im.ln(); }, || { let _ = re_im; }, batch, reps)));
    /// ```
    #[inline(never)]
    pub fn ln(self) -> Self {
        use Cpx::*;
        match self {
            One => Zero,
            NegOne => Im { im: T::PI },
            J => Im { im: T::FRAC_PI_2 },
            NegJ => Im {
                im: T::NEG_FRAC_PI_2,
            },
            Zero => Re {
                re: T::NEG_INFINITY,
            },
            Ph { ph: im } => Im { im },
            Ln { re, im } => ReIm { re, im },
            Pl { rad, ph: im } => ReIm { re: rad.ln(), im },
            Re { re } => match re {
                x if x > T::ZERO => Re { re: x.ln() },
                x if x < T::ZERO => ReIm {
                    re: (-x).ln(),
                    im: T::PI,
                },
                _ => Zero,
            },
            Im { im } => match im {
                x if x > T::ZERO => ReIm {
                    re: x.ln(),
                    im: T::FRAC_PI_2,
                },
                x if x < T::ZERO => ReIm {
                    re: (-x).ln(),
                    im: T::NEG_FRAC_PI_2,
                },
                _ => Zero,
            },
            ReIm { re, im } => ReIm {
                re: re.hypot(im).ln(),
                im: im.atan2(re),
            },
        }
    }
    /// Raises the complex number to an integer power using polar form.
    ///
    /// Computes `rad^n * e^(i * n * θ)` by raising the magnitude to `n` and scaling the phase.
    /// Returns the canonicalized result.
    ///
    /// # Example
    /// ```rust
    /// use cpx_coords::{Cpx,FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (v0, v1): (F,F) = (3.0,4.0);
    ///
    /// // Any number to the 0th power is 1
    /// let z: C = ReIm { re: v0, im: v1 };
    /// //assert_eq!(z.powi32(0), One);
    ///
    /// let pow: i32 = 2;
    /// // Real positive number
    /// let x: C = Re { re: v0 };
    /// //assert_eq!(x.powi32(pow), Re {re: v0.powi(pow)});
    ///
    /// // Real negative number to even power: positive real
    /// let x: C = Re { re: -v0 };
    /// //assert_eq!(x.powi32(pow), Re {re: (-v0).powi(pow)});
    ///
    /// // Imaginary unit squared: i^2 = -1
    /// //assert_eq!(C::J.powi32(2), NegOne);
    ///
    /// // (i)^4 = 1
    /// //assert_eq!(C::J.powi32(4), One);
    ///
    /// let phase: F = F::FRAC_PI_2;
    ///
    /// // Polar representation
    /// let z: Cpx<F> = Pl { rad: v0, ph: phase };
    /// //assert_eq!(z.powi32(pow), Pl { rad: v0.powi(pow), ph: phase * F::from(pow)});
    ///
    /// // General complex number: verify magnitude grows appropriately
    /// let z: C = ReIm { re: 1.0, im: 1.0 };
    /// let z2 = z.powi32(2);
    /// //assert_eq!(z.powi32(2), Im { im: 2.0}); // This fails due to floating-point precision limits in f32/f64.
    /// //assert!((z2 - Im { im: 2.0}).rad() < F::THRESHOLD);
    ///
    /// // Zero to positive power is still T::ZERO
    /// //assert_eq!(C::Zero.powi32(5), Zero);
    ///
    /// // Zero to T::ZERO is defined as One (by convention)
    /// //assert_eq!(C::Zero.powi32(0), One);
    ///
    /// // J = i; i^3 = -i = -J
    /// //assert_eq!(C::J.powi32(3), NegJ);
    ///
    /// // Negative J to odd power remains imaginary
    /// //assert_eq!(C::NegJ.powi32(3), J);
    ///
    /// // Unit modulus phase rotation (θ = π/2): (e^{iπ/2})^4 = 1
    /// let z = Ph { ph: F::FRAC_PI_2 };
    /// //assert_eq!(z.powi32(4), One);
    ///
    /// // Phase angle wrapping: e^{iπ} = -1; (e^{iπ})^3 = -1
    /// let z = Ph { ph: F::PI };
    /// //assert_eq!(z.powi32(3), NegOne);
    ///
    /// // Large integer exponent: (2 + 0i)^20 = 2^20
    /// let x: C = Re { re: 2.0 };
    /// //assert_eq!(x.powi32(20), Re { re: 1048576.0 });
    /// ```
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 10;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = zero.powi32(2); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = one.powi32(2); }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=7).contains(&measure_cycles( || { let _ = neg_one.powi32(2); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=8).contains(&measure_cycles( || { let _ = j.powi32(2); }, || { let _ = j; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = neg_j.powi32(2); }, || { let _ = neg_j; }, batch, reps)));
    /// // assert!((..=25).contains(&measure_cycles( || { let _ = re.powi32(2); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=30).contains(&measure_cycles( || { let _ = im.powi32(2); }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=28).contains(&measure_cycles( || { let _ = ph.powi32(2); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = pl.powi32(2); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = ln.powi32(2); }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = re_im.powi32(0); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = re_im.powi32(1); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=45).contains(&measure_cycles( || { let _ = re_im.powi32(-1); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=30).contains(&measure_cycles( || { let _ = re_im.powi32(2); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=50).contains(&measure_cycles( || { let _ = re_im.powi32(-2); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=55).contains(&measure_cycles( || { let _ = re_im.powi32(3); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=50).contains(&measure_cycles( || { let _ = re_im.powi32(4); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=70).contains(&measure_cycles( || { let _ = re_im.powi32(8); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=90).contains(&measure_cycles( || { let _ = re_im.powi32(16); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=110).contains(&measure_cycles( || { let _ = re_im.powi32(32); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=130).contains(&measure_cycles( || { let _ = re_im.powi32(64); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=75).contains(&measure_cycles( || { let _ = re_im.powi32(6); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=95).contains(&measure_cycles( || { let _ = re_im.powi32(10); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=95).contains(&measure_cycles( || { let _ = re_im.powi32(12); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=115).contains(&measure_cycles( || { let _ = re_im.powi32(14); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=115).contains(&measure_cycles( || { let _ = re_im.powi32(18); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=115).contains(&measure_cycles( || { let _ = re_im.powi32(20); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=130).contains(&measure_cycles( || { let _ = re_im.powi32(22); }, || { let _ = re_im; }, batch, reps))); // No benefit
    /// // assert!((..=115).contains(&measure_cycles( || { let _ = re_im.powi32(24); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=155).contains(&measure_cycles( || { let _ = re_im.powi32(1000); }, || { let _ = re_im; }, batch, reps)));
    /// ```
    #[inline(never)]
    pub fn powi32(self, n: i32) -> Self {
        use Cpx::*;

        // Special cases first
        if n == 0 {
            return One;
        } else if n == 1 {
            return self;
        } else if n == -1 {
            return self.recip();
        }
        #[inline(always)]
        fn zero_to_pow<T: FloatExt>(n: i32) -> Cpx<T> {
            match n {
                _ if n > 0 => Zero,
                _ => Re { re: T::INFINITY },
            }
        }

        #[inline(always)]
        fn neg_one_to_pow<T: FloatExt>(n: i32) -> Cpx<T> {
            match n & 1_i32 {
                0 => One,
                _ => NegOne,
            }
        }

        #[inline(always)]
        fn j_to_pow<T: FloatExt>(n: i32) -> Cpx<T> {
            match n & 3_i32 {
                0 => One,
                1 => J,
                2 => NegOne,
                _ => NegJ,
            }
        }

        #[inline(always)]
        fn neg_j_to_pow<T: FloatExt>(n: i32) -> Cpx<T> {
            match n & 3_i32 {
                0 => One,
                1 => NegJ,
                2 => NegOne,
                _ => J,
            }
        }

        #[inline(always)]
        fn ph_times_i32<T: FloatExt>(ph: T, n: i32) -> T {
            (ph * T::from(n).unwrap()).principal_val()
        }

        match self {
            One => One,
            Zero => zero_to_pow(n),
            NegOne => neg_one_to_pow(n),
            J => j_to_pow(n),
            NegJ => neg_j_to_pow(n),
            Re { re } => Re { re: re.powi(n) },
            Im { im } => {
                let val = im.powi(n);
                match n & 3_i32 {
                    0 => Re { re: val },
                    1 => Im { im: val },
                    2 => Re { re: -val },
                    _ => Im { im: -val },
                }
            }
            Ph { ph } => Ph {
                ph: ph_times_i32(ph, n),
            },
            Pl { rad, ph } => Pl {
                rad: rad.powi(n),
                ph: ph_times_i32(ph, n),
            },
            Ln { re, im: ph } => Ln {
                re: re * T::from(n).unwrap(),
                im: ph_times_i32(ph, n),
            },
            ReIm { re, im } => {
                #[inline(always)]
                fn square<T: FloatExt>(re: T, im: T) -> (T, T) {
                    let re2 = re * re - im * im;
                    let im2 = re * im * T::TWO;
                    (re2, im2)
                }

                #[inline(always)]
                fn tiny_mul<T: FloatExt>(lhs: (T, T), rhs: (T, T)) -> (T, T) {
                    let re = lhs.0 * rhs.0 - lhs.1 * rhs.1;
                    let im = lhs.1 * rhs.0 + lhs.0 * rhs.1;
                    (re, im)
                }

                match n {
                    2 => {
                        let (re_2, im_2) = square(re, im);
                        ReIm { re: re_2, im: im_2 }
                    }
                    4 => {
                        let (re_2, im_2) = square(re, im); // z^2
                        let (re_4, im_4) = square(re_2, im_2); // (z^2)^2 = z^4
                        ReIm { re: re_4, im: im_4 }
                    }
                    8 => {
                        let (re_2, im_2) = square(re, im); // z^2
                        let (re_4, im_4) = square(re_2, im_2); // (z^2)^2 = z^4
                        let (re_8, im_8) = square(re_4, im_4); // (z^4)^2 = z^8
                        ReIm { re: re_8, im: im_8 }
                    }
                    16 => {
                        let (re_2, im_2) = square(re, im);
                        let (re_4, im_4) = square(re_2, im_2);
                        let (re_8, im_8) = square(re_4, im_4);
                        let (re_16, im_16) = square(re_8, im_8);
                        ReIm {
                            re: re_16,
                            im: im_16,
                        }
                    }
                    32 => {
                        let (re_2, im_2) = square(re, im);
                        let (re_4, im_4) = square(re_2, im_2);
                        let (re_8, im_8) = square(re_4, im_4);
                        let (re_16, im_16) = square(re_8, im_8);
                        let (re_32, im_32) = square(re_16, im_16);
                        ReIm {
                            re: re_32,
                            im: im_32,
                        }
                    }
                    -2 => {
                        let (re_2, im_2) = square(re, im);
                        let rad2 = re * re + im * im;
                        let rad4 = rad2 * rad2;
                        ReIm {
                            re: re_2 / rad4,
                            im: -im_2 / rad4,
                        }
                    }
                    3 => {
                        let (re_2, im_2) = square(re, im);
                        let (re, im) = tiny_mul((re_2, im_2), (re, im));
                        ReIm { re, im }
                    }
                    6 => {
                        let (re_2, im_2) = square(re, im);
                        let (re_4, im_4) = square(re_2, im_2);
                        let (re, im) = tiny_mul((re_4, im_4), (re_2, im_2));
                        ReIm { re, im }
                    }
                    10 => {
                        let (re_2, im_2) = square(re, im);
                        let (re_4, im_4) = square(re_2, im_2);
                        let (re_8, im_8) = square(re_4, im_4);
                        let (re, im) = tiny_mul((re_8, im_8), (re_2, im_2));
                        ReIm { re, im }
                    }
                    12 => {
                        let (re_2, im_2) = square(re, im);
                        let (re_4, im_4) = square(re_2, im_2);
                        let (re_8, im_8) = square(re_4, im_4);
                        let (re, im) = tiny_mul((re_8, im_8), (re_4, im_4));
                        ReIm { re, im }
                    }
                    14 => {
                        let (re_2, im_2) = square(re, im);
                        let (re_4, im_4) = square(re_2, im_2);
                        let (re_8, im_8) = square(re_4, im_4);
                        let (re_12, im_12) = tiny_mul((re_8, im_8), (re_4, im_4));
                        let (re, im) = tiny_mul((re_12, im_12), (re_2, im_2));
                        ReIm { re, im }
                    }
                    18 => {
                        let (re_2, im_2) = square(re, im);
                        let (re_4, im_4) = square(re_2, im_2);
                        let (re_8, im_8) = square(re_4, im_4);
                        let (re_16, im_16) = square(re_8, im_8);
                        let (re, im) = tiny_mul((re_16, im_16), (re_2, im_2));
                        ReIm { re, im }
                    }
                    20 => {
                        let (re_2, im_2) = square(re, im);
                        let (re_4, im_4) = square(re_2, im_2);
                        let (re_8, im_8) = square(re_4, im_4);
                        let (re_16, im_16) = square(re_8, im_8);
                        let (re, im) = tiny_mul((re_16, im_16), (re_4, im_4));
                        ReIm { re, im }
                    }
                    //22 => {
                    //    let (re_2, im_2) = square(re, im);
                    //    let (re_4, im_4) = square(re_2, im_2);
                    //    let (re_8, im_8) = square(re_4, im_4);
                    //    let (re_16, im_16) = square(re_8, im_8);
                    //    let (re_20, im_20) = tiny_mul((re_16, im_16), (re_4, im_4));
                    //    let (re, im) = tiny_mul((re_20, im_20), (re_2, im_2));
                    //    ReIm { re, im }
                    //}
                    24 => {
                        let (re_2, im_2) = square(re, im);
                        let (re_4, im_4) = square(re_2, im_2);
                        let (re_8, im_8) = square(re_4, im_4);
                        let (re_16, im_16) = square(re_8, im_8);
                        let (re, im) = tiny_mul((re_16, im_16), (re_8, im_8));
                        ReIm { re, im }
                    }
                    _ => {
                        let (rad, ph) = (re.hypot(im), im.atan2(re)); // self.rad_ph();
                        Pl {
                            rad: rad.powi(n),
                            ph: ph_times_i32(ph, n),
                        }
                    }
                }
            }
        }
    }

    /// Raises the complex number to a floating-point power using polar form.
    ///
    /// Computes `r^f * e^(i * f * θ)` by raising the magnitude to `f` and scaling the phase.
    /// Returns the canonicalized result.
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=11).contains(&measure_cycles( || { let _ = zero.powf(0.0); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=25).contains(&measure_cycles( || { let _ = one.powf(0.5); }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=25).contains(&measure_cycles( || { let _ = zero.powf(0.5); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = neg_one.powf(0.5); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = j.powf(0.5); }, || { let _ = j; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = neg_j.powf(0.5); }, || { let _ = neg_j; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = ph.powf(0.5); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = ln.powf(0.5); }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=50).contains(&measure_cycles( || { let _ = re.powf(0.5); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=60).contains(&measure_cycles( || { let _ = pl.powf(0.5); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=60).contains(&measure_cycles( || { let _ = im.powf(0.5); }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=110).contains(&measure_cycles( || { let _ = re_im.powf(0.5); }, || { let _ = re_im; }, batch, reps)));
    /// ```
    pub fn powf(self, f: T) -> Self {
        use Cpx::*;

        if f == T::ZERO {
            return One;
        } else if f == T::ONE {
            return self;
        } else if f == T::NEG_ONE {
            return self.recip();
        }
        #[inline(always)]
        fn ph_times_f<T: FloatExt>(ph: T, f: T) -> T {
            (ph * f).principal_val()
        }

        #[inline(always)]
        fn f_mod2<T: FloatExt>(f: T) -> T {
            if f >= -T::TWO && f <= T::TWO {
                f
            } else {
                f % T::TWO
            }
        }

        #[inline(always)]
        fn f_mod4<T: FloatExt>(f: T) -> T {
            let four = T::TWO + T::TWO;
            if f >= -four && f <= four { f } else { f % four }
        }

        match self {
            One => One,
            Zero => match f {
                _ if f > T::ZERO => Zero,
                _ => Re { re: T::INFINITY },
            },
            NegOne => Ph {
                ph: ph_times_f(T::PI, f_mod2(f)),
            },
            J => Ph {
                ph: ph_times_f(T::FRAC_PI_2, f_mod4(f)),
            },
            NegJ => Ph {
                ph: ph_times_f(T::NEG_FRAC_PI_2, f_mod4(f)),
            },
            Ph { ph } => Ph {
                ph: ph_times_f(ph, f),
            },
            Ln { re, im: ph } => Ln {
                re: re * f,
                im: ph_times_f(ph, f),
            },
            Pl { rad, ph } => Pl {
                rad: rad.powf(f),
                ph: ph_times_f(ph, f),
            },
            Re { re } => match re {
                _ if re > T::ZERO => Re { re: re.powf(f) },
                _ => Pl {
                    rad: (-re).powf(f),
                    ph: ph_times_f(T::PI, f_mod2(f)),
                },
            },
            Im { im } => match im {
                _ if im > T::ZERO => Pl {
                    rad: im.powf(f),
                    ph: ph_times_f(T::FRAC_PI_2, f_mod4(f)),
                },
                _ => Pl {
                    rad: (-im).powf(f),
                    ph: ph_times_f(T::NEG_FRAC_PI_2, f_mod4(f)),
                },
            },
            ReIm { re, im } => Pl {
                rad: (re * re + im * im).powf(f / T::TWO),
                ph: ph_times_f(im.atan2(re), f),
            },
        }
    }
    /// Returns the *n*-th root of the complex number.
    ///
    /// This computes `self^(1/n)` using floating-point exponentiation in polar form.
    ///
    /// # Panics
    ///
    /// Panics if `n == 0` (root index cannot be zero) or if the conversion from `u32` to `T` fails.
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = zero.root_u32(1); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=20).contains(&measure_cycles( || { let _ = zero.root_u32(2); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=20).contains(&measure_cycles( || { let _ = one.root_u32(2); }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=90).contains(&measure_cycles( || { let _ = neg_one.root_u32(2); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=90).contains(&measure_cycles( || { let _ = j.root_u32(2); }, || { let _ = j; }, batch, reps)));
    /// // assert!((..=90).contains(&measure_cycles( || { let _ = neg_j.root_u32(2); }, || { let _ = neg_j; }, batch, reps)));
    /// // assert!((..=80).contains(&measure_cycles( || { let _ = ph.root_u32(2); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=190).contains(&measure_cycles( || { let _ = ln.root_u32(2); }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=130).contains(&measure_cycles( || { let _ = re.root_u32(2); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=140).contains(&measure_cycles( || { let _ = pl.root_u32(2); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=140).contains(&measure_cycles( || { let _ = im.root_u32(2); }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=170).contains(&measure_cycles( || { let _ = re_im.root_u32(2); }, || { let _ = re_im; }, batch, reps)));
    /// ```
    pub fn root_u32(self, n: u32) -> Self {
        use Cpx::*;
        if n == 1_u32 {
            return self;
        } else if n == 0_u32 {
            panic!("root number cannot be zero");
        };

        let canon = self.canonicalize_lazy();

        match canon {
            Zero | One => canon,
            _ => {
                // Precompute conversions once
                let n_t: T = T::from(n).unwrap();
                let recip = n_t.recip();
                canon.powf(recip)
            }
        }
    }
    /// Raises the complex number to a complex exponent `cpx` using the polar form.
    ///
    /// This computes
    ///
    /// ```text
    /// z^c = r^c * e^{i * c * θ}
    /// ```
    ///
    /// where `r` and `θ` are the magnitude and phase of `z`.
    /// The result is returned in canonical form.
    ///
    /// # Mathematical note
    ///
    /// If `z = 0` and `c = a + i b`:
    /// - If `b ≠ 0`, the limit `lim_{z→0} z^c` does not exist;
    /// - If `b = 0` and `a > 0`, the limit is `0`;
    /// - If `b = 0` and `a = 0`, the limit is `1`;
    /// - If `b = 0` and `a < 0`, the limit diverges to `+∞`.
    ///
    /// This implementation panics for the non-existent case (`b ≠ 0`) and for the divergent case (`a < 0`).
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (mut zero, mut one, mut neg_one, mut j, mut neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (mut re, mut im, mut re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (mut ph, mut pl, mut ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = re_im.powc(zero); }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=8).contains(&measure_cycles( || { let _ = one.powc(re_im); }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=20).contains(&measure_cycles( || { let _ = zero.powc(one); }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=15).contains(&measure_cycles( || { let _ = neg_one.powc(neg_one); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=28).contains(&measure_cycles( || { let _ = neg_one.powc(j); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=30).contains(&measure_cycles( || { let _ = neg_one.powc(re); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=65).contains(&measure_cycles( || { let _ = neg_one.powc(im); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=35).contains(&measure_cycles( || { let _ = neg_one.powc(re_im); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=110).contains(&measure_cycles( || { let _ = neg_one.powc(ph); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=90).contains(&measure_cycles( || { let _ = neg_one.powc(pl); }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=30).contains(&measure_cycles( || { let _ = j.powc(re); }, || { let _ = j; }, batch, reps)));
    /// // assert!((..=65).contains(&measure_cycles( || { let _ = neg_j.powc(im); }, || { let _ = neg_j; }, batch, reps)));
    /// // assert!((..=35).contains(&measure_cycles( || { let _ = ph.powc(re_im); }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=45).contains(&measure_cycles( || { let _ = ln.powc(re_im); }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=45).contains(&measure_cycles( || { let _ = re.powc(re_im); }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=65).contains(&measure_cycles( || { let _ = pl.powc(re_im); }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=60).contains(&measure_cycles( || { let _ = im.powc(re_im); }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=150).contains(&measure_cycles( || { let _ = re_im.powc(re_im); }, || { let _ = re_im; }, batch, reps)));
    /// ```
    pub fn powc(self, c: Self) -> Self {
        use Cpx::*;
        match (self, c) {
            (_, Zero) | (One, _) => One,
            (_, One) => self,
            (_, NegOne) => self.recip(),
            (Zero, x) => match x {
                Zero => One,
                One => Zero,
                NegOne => Re { re: T::INFINITY },
                Re { re } => match re {
                    x if x > T::ZERO => Zero,
                    x if x < T::ZERO => Re { re: T::INFINITY },
                    _ => One,
                },
                Im { im } => {
                    if im == T::ZERO {
                        One
                    } else {
                        panic!("Zero to the complex power is undefined")
                    }
                }
                J | NegJ => panic!("Zero to the complex power is undefined"),
                _ => {
                    let (re, im) = c.re_im();
                    match im {
                        x if x == T::ZERO => match re {
                            x if x > T::ZERO => Zero,
                            x if x < T::ZERO => Re { re: T::INFINITY },
                            _ => One,
                        },
                        _ => panic!("Zero to the complex power is undefined"),
                    }
                }
            },
            (NegOne, x) => match x {
                J => Re { re: T::EXP_NEG_PI },
                NegJ => Re { re: T::EXP_PI },
                Re { re } => Ph { ph: T::PI * re },
                Im { im } => Re {
                    re: T::EXP_NEG_PI.powf(im),
                },
                ReIm { re, im } => Ln {
                    re: -T::PI * im,
                    im: T::PI * re,
                },
                Ph { ph } => Ln {
                    re: -T::PI * ph.sin(),
                    im: T::PI * ph.cos(),
                },
                Pl { rad, ph } => {
                    let new_rad = T::PI * rad;
                    Ln {
                        re: -new_rad * ph.sin(),
                        im: new_rad * ph.cos(),
                    }
                }
                Ln { re, im } => {
                    let rad = T::PI * re.exp();
                    Ln {
                        re: -rad * im.sin(),
                        im: rad * im.cos(),
                    }
                }
                _ => unreachable!(),
            },
            (J, x) => match x {
                J => Re {
                    re: T::EXP_NEG_FRAC_PI_2,
                },
                NegJ => Re {
                    re: T::EXP_FRAC_PI_2,
                },
                Re { re } => Ph {
                    ph: T::FRAC_PI_2 * re,
                },
                Im { im } => Re {
                    re: T::EXP_NEG_FRAC_PI_2.powf(im),
                },
                ReIm { re, im } => Ln {
                    re: T::NEG_FRAC_PI_2 * im,
                    im: T::FRAC_PI_2 * re,
                },
                Ph { ph } => Ln {
                    re: T::NEG_FRAC_PI_2 * ph.sin(),
                    im: T::FRAC_PI_2 * ph.cos(),
                },
                Pl { rad, ph } => {
                    let (r2, i2) = (rad * ph.cos(), rad * ph.sin());
                    Ln {
                        re: T::NEG_FRAC_PI_2 * i2,
                        im: T::FRAC_PI_2 * r2,
                    }
                }
                Ln { re, im } => {
                    let rad = re.exp();
                    let (r2, i2) = (rad * im.cos(), rad * im.sin());
                    Ln {
                        re: T::NEG_FRAC_PI_2 * i2,
                        im: T::FRAC_PI_2 * r2,
                    }
                }
                _ => unreachable!(),
            },
            (NegJ, x) => match x {
                J => Re {
                    re: T::EXP_FRAC_PI_2,
                },
                NegJ => Re {
                    re: T::EXP_NEG_FRAC_PI_2,
                },
                Re { re } => Ph {
                    ph: T::NEG_FRAC_PI_2 * re,
                },
                Im { im } => Re {
                    re: T::EXP_FRAC_PI_2.powf(im),
                },
                ReIm { re, im } => Ln {
                    re: T::FRAC_PI_2 * im,
                    im: T::NEG_FRAC_PI_2 * re,
                },
                Ph { ph } => Ln {
                    re: T::FRAC_PI_2 * ph.sin(),
                    im: T::NEG_FRAC_PI_2 * ph.cos(),
                },
                Pl { rad, ph } => {
                    let (r2, i2) = (rad * ph.cos(), rad * ph.sin());
                    Ln {
                        re: T::FRAC_PI_2 * i2,
                        im: T::NEG_FRAC_PI_2 * r2,
                    }
                }
                Ln { re, im } => {
                    let rad = re.exp();
                    let (r2, i2) = (rad * im.cos(), rad * im.sin());
                    Ln {
                        re: T::FRAC_PI_2 * i2,
                        im: T::NEG_FRAC_PI_2 * r2,
                    }
                }
                _ => unreachable!(),
            },
            (Ph { ph }, x) => {
                let new_ph = ph.principal_val();
                match x {
                    J => Re {
                        re: (-new_ph).exp(),
                    },
                    NegJ => Re { re: new_ph.exp() },
                    Re { re } => Ph { ph: new_ph * re },
                    Im { im } => Re {
                        re: (-new_ph * im).exp(),
                    },
                    ReIm { re, im } => Ln {
                        re: -new_ph * im,
                        im: new_ph * re,
                    },
                    Ph { ph: ph2 } => Ln {
                        re: -new_ph * ph2.sin(),
                        im: new_ph * ph2.cos(),
                    },
                    Pl { rad, ph: ph2 } => {
                        let (r2, i2) = (rad * ph2.cos(), rad * ph2.sin());
                        Ln {
                            re: -new_ph * i2,
                            im: new_ph * r2,
                        }
                    }
                    Ln { re, im } => {
                        let rad = re.exp();
                        let (r2, i2) = (rad * im.cos(), rad * im.sin());
                        Ln {
                            re: -new_ph * i2,
                            im: new_ph * r2,
                        }
                    }
                    _ => unreachable!(),
                }
            }
            (Ln { re, im }, x) => {
                let im = im.principal_val();
                match x {
                    J => Ln { re: -im, im: re },
                    NegJ => Ln { re: im, im: -re },
                    Re { re: r2 } => Ln {
                        re: re * r2,
                        im: im * r2,
                    },
                    Im { im: i2 } => Ln {
                        re: -im * i2,
                        im: re * i2,
                    },
                    ReIm { re: r2, im: i2 } => Ln {
                        re: re * r2 - im * i2,
                        im: re * i2 + im * r2,
                    },
                    Ph { ph } => {
                        let (r2, i2) = (ph.cos(), ph.sin());
                        Ln {
                            re: re * r2 - im * i2,
                            im: re * i2 + im * r2,
                        }
                    }
                    Pl { rad, ph } => {
                        let (r2, i2) = (rad * ph.cos(), rad * ph.sin());
                        Ln {
                            re: re * r2 - im * i2,
                            im: re * i2 + im * r2,
                        }
                    }
                    Ln { re: r3, im: ph } => {
                        let rad = r3.exp();
                        let (r2, i2) = (rad * ph.cos(), rad * ph.sin());
                        Ln {
                            re: re * r2 - im * i2,
                            im: re * i2 + im * r2,
                        }
                    }
                    _ => unreachable!(),
                }
            }
            (Pl { rad, ph }, x) => {
                let ph = ph.principal_val();
                let ln_rad = rad.ln();
                match x {
                    J => Ln {
                        re: -ph,
                        im: ln_rad,
                    },
                    NegJ => Ln {
                        re: ph,
                        im: -ln_rad,
                    },
                    Re { re } => Ln {
                        re: ln_rad * re,
                        im: ph * re,
                    },
                    Im { im } => Ln {
                        re: -ph * im,
                        im: ln_rad * im,
                    },
                    ReIm { re, im } => Ln {
                        re: ln_rad * re - ph * im,
                        im: ln_rad * im + ph * re,
                    },
                    Ph { ph: ph2 } => {
                        let (r2, i2) = (ph2.cos(), ph2.sin());
                        Ln {
                            re: ln_rad * r2 - ph * i2,
                            im: ln_rad * i2 + ph * r2,
                        }
                    }
                    Pl { rad: rad2, ph: ph2 } => {
                        let (r2, i2) = (rad2 * ph2.cos(), rad2 * ph2.sin());
                        Ln {
                            re: ln_rad * r2 - ph * i2,
                            im: ln_rad * i2 + ph * r2,
                        }
                    }
                    Ln { re, im } => {
                        let rad = re.exp();
                        let (r2, i2) = (rad * im.cos(), rad * im.sin());
                        Ln {
                            re: ln_rad * r2 - ph * i2,
                            im: ln_rad * i2 + ph * r2,
                        }
                    }
                    _ => unreachable!(),
                }
            }
            (Re { re }, x) => match re {
                _ if re > T::ZERO => {
                    let ln_rad = re.ln();
                    match x {
                        J => Ph { ph: ln_rad },
                        NegJ => Ph { ph: -ln_rad },
                        Re { re: r2 } => Re { re: re.powf(r2) },
                        Im { im } => Ph { ph: ln_rad * im },
                        ReIm { re: r2, im } => Ln {
                            re: ln_rad * r2,
                            im: ln_rad * im,
                        },
                        Ph { ph } => {
                            let (r2, i2) = (ph.cos(), ph.sin());
                            Ln {
                                re: ln_rad * r2,
                                im: ln_rad * i2,
                            }
                        }
                        Pl { rad, ph } => {
                            let (r2, i2) = (rad * ph.cos(), rad * ph.sin());
                            Ln {
                                re: ln_rad * r2,
                                im: ln_rad * i2,
                            }
                        }
                        Ln { re: r3, im } => {
                            let rad = r3.exp();
                            let (r2, i2) = (rad * im.cos(), rad * im.sin());
                            Ln {
                                re: ln_rad * r2,
                                im: ln_rad * i2,
                            }
                        }
                        _ => unreachable!(),
                    }
                }
                _ if re < T::ZERO => {
                    let ln_rad = (-re).ln();
                    match x {
                        J => Pl {
                            rad: T::EXP_NEG_PI,
                            ph: ln_rad,
                        },
                        NegJ => Pl {
                            rad: T::EXP_PI,
                            ph: -ln_rad,
                        },
                        Re { re: r2 } => Ln {
                            re: ln_rad * r2,
                            im: T::PI * r2,
                        },
                        Im { im } => Ln {
                            re: -T::PI * im,
                            im: ln_rad * im,
                        },
                        ReIm { re: r2, im } => Ln {
                            re: ln_rad * r2 + -T::PI * im,
                            im: T::PI * r2 + ln_rad * im,
                        },
                        Ph { ph } => {
                            let (r2, im) = (ph.cos(), ph.sin());
                            Ln {
                                re: ln_rad * r2 + -T::PI * im,
                                im: T::PI * r2 + ln_rad * im,
                            }
                        }
                        Pl { rad, ph } => {
                            let (r2, im) = (rad * ph.cos(), rad * ph.sin());
                            Ln {
                                re: ln_rad * r2 + -T::PI * im,
                                im: T::PI * r2 + ln_rad * im,
                            }
                        }
                        Ln { re: r3, im } => {
                            let rad = r3.exp();
                            let (r2, im) = (rad * im.cos(), rad * im.sin());
                            Ln {
                                re: ln_rad * r2 + -T::PI * im,
                                im: T::PI * r2 + ln_rad * im,
                            }
                        }
                        _ => unreachable!(),
                    }
                }
                _ => match x {
                    J | NegJ => panic!("Zero to the complex power is undefined"),
                    Re { re: r2 } => match r2 {
                        _ if r2 > T::ZERO => Zero,
                        _ if r2 < T::ZERO => Re { re: T::infinity() },
                        _ => One,
                    },
                    Im { im } => match im {
                        _ if im == T::ZERO => One,
                        _ => panic!("Zero to the complex power is undefined"),
                    },
                    ReIm { re: r2, im } => match im {
                        _ if im == T::ZERO => match r2 {
                            _ if r2 > T::ZERO => Zero,
                            _ if r2 < T::ZERO => Re { re: T::infinity() },
                            _ => One,
                        },
                        _ => panic!("Zero to the complex power is undefined"),
                    },
                    Ph { ph } => match ph.principal_val() {
                        x if x == T::ZERO => Zero,
                        x if x == T::PI => Re { re: T::infinity() },
                        _ => panic!("Zero to the complex power is undefined"),
                    },
                    Pl { rad, ph } => match rad {
                        _ if rad == T::ZERO => One,
                        _ => match ph.principal_val() {
                            x if x == T::ZERO => Zero,
                            x if x == T::PI => Re { re: T::infinity() },
                            _ => panic!("Zero to the complex power is undefined"),
                        },
                    },
                    Ln { re: r2, im: ph } => match r2 {
                        _ if r2 < T::LN_THRESHOLD => One,
                        _ => match ph.principal_val() {
                            x if x == T::ZERO => Zero,
                            x if x == T::PI => Re { re: T::infinity() },
                            _ => panic!("Zero to the complex power is undefined"),
                        },
                    },
                    _ => unreachable!(),
                },
            },
            (Im { im }, x) => match im {
                _ if im > T::ZERO => {
                    let ln_im = im.ln();
                    match x {
                        J => Pl {
                            rad: T::EXP_NEG_FRAC_PI_2,
                            ph: ln_im,
                        },
                        NegJ => Pl {
                            rad: T::EXP_FRAC_PI_2,
                            ph: -ln_im,
                        },
                        Re { re } => Ln {
                            re: ln_im * re,
                            im: T::FRAC_PI_2 * re,
                        },
                        Im { im: i2 } => Ln {
                            re: T::NEG_FRAC_PI_2 * i2,
                            im: ln_im * i2,
                        },
                        ReIm { re, im: i2 } => Ln {
                            re: ln_im * re + T::NEG_FRAC_PI_2 * i2,
                            im: ln_im * i2 + T::FRAC_PI_2 * re,
                        },
                        Ph { ph } => {
                            let (re, i2) = (ph.cos(), ph.sin());
                            Ln {
                                re: ln_im * re + T::NEG_FRAC_PI_2 * i2,
                                im: ln_im * i2 + T::FRAC_PI_2 * re,
                            }
                        }
                        Pl { rad, ph } => {
                            let (re, i2) = (rad * ph.cos(), rad * ph.sin());
                            Ln {
                                re: ln_im * re + T::NEG_FRAC_PI_2 * i2,
                                im: ln_im * i2 + T::FRAC_PI_2 * re,
                            }
                        }
                        Ln { re, im: i3 } => {
                            let rad = re.exp();
                            let (re, i2) = (rad * i3.cos(), rad * i3.sin());
                            Ln {
                                re: ln_im * re + T::NEG_FRAC_PI_2 * i2,
                                im: ln_im * i2 + T::FRAC_PI_2 * re,
                            }
                        }
                        _ => unreachable!(),
                    }
                }
                _ if im < T::ZERO => {
                    let ln_im = (-im).ln();
                    match x {
                        J => Pl {
                            rad: T::EXP_FRAC_PI_2,
                            ph: -ln_im,
                        },
                        NegJ => Pl {
                            rad: T::EXP_NEG_FRAC_PI_2,
                            ph: ln_im,
                        },
                        Re { re } => Ln {
                            re: ln_im * re,
                            im: T::NEG_FRAC_PI_2 * re,
                        },
                        Im { im: i2 } => Ln {
                            re: T::FRAC_PI_2 * i2,
                            im: ln_im * i2,
                        },
                        ReIm { re, im: i2 } => Ln {
                            re: ln_im * re + T::FRAC_PI_2 * i2,
                            im: ln_im * i2 + T::NEG_FRAC_PI_2 * re,
                        },
                        Ph { ph } => {
                            let (re, i2) = (ph.cos(), ph.sin());
                            Ln {
                                re: ln_im * re + T::FRAC_PI_2 * i2,
                                im: ln_im * i2 + T::NEG_FRAC_PI_2 * re,
                            }
                        }
                        Pl { rad, ph } => {
                            let (re, i2) = (rad * ph.cos(), rad * ph.sin());
                            Ln {
                                re: ln_im * re + T::FRAC_PI_2 * i2,
                                im: ln_im * i2 + T::NEG_FRAC_PI_2 * re,
                            }
                        }
                        Ln { re, im: ph } => {
                            let rad = re.exp();
                            let (re, i2) = (rad * ph.cos(), rad * ph.sin());
                            Ln {
                                re: ln_im * re + T::FRAC_PI_2 * i2,
                                im: ln_im * i2 + T::NEG_FRAC_PI_2 * re,
                            }
                        }
                        _ => unreachable!(),
                    }
                }
                _ => match x {
                    J | NegJ => panic!("Zero to the complex power is undefined"),
                    Re { re } => match re {
                        _ if re > T::ZERO => Zero,
                        _ if re < T::ZERO => Re { re: T::infinity() },
                        _ => One,
                    },
                    Im { im: i2 } => match i2 {
                        _ if im == T::ZERO => One,
                        _ => panic!("Zero to the complex power is undefined"),
                    },
                    ReIm { re, im: i2 } => match i2 {
                        _ if i2 == T::ZERO => match re {
                            _ if re > T::ZERO => Zero,
                            _ if re < T::ZERO => Re { re: T::infinity() },
                            _ => One,
                        },
                        _ => panic!("Zero to the complex power is undefined"),
                    },
                    Ph { ph } => match ph.principal_val() {
                        x if x == T::ZERO => Zero,
                        x if x == T::PI => Re { re: T::infinity() },
                        _ => panic!("Zero to the complex power is undefined"),
                    },
                    Pl { rad, ph } => match rad {
                        x if x == T::ZERO => One,
                        _ => match ph.principal_val() {
                            x if x == T::ZERO => Zero,
                            x if x == T::PI => Re { re: T::infinity() },
                            _ => panic!("Zero to the complex power is undefined"),
                        },
                    },
                    Ln { re, im: i2 } => match re {
                        x if x < T::LN_THRESHOLD => One,
                        _ => match i2.principal_val() {
                            x if x == T::ZERO => Zero,
                            x if x == T::PI => Re { re: T::infinity() },
                            _ => panic!("Zero to the complex power is undefined"),
                        },
                    },
                    _ => unreachable!(),
                },
            },
            (ReIm { re, im }, x) => {
                let rad = re.hypot(im);
                match rad {
                    _ if rad == T::ZERO => match x {
                        J | NegJ => panic!("Zero to the complex power is undefined"),
                        Re { re } => match re {
                            _ if re > T::ZERO => Zero,
                            _ if re < T::ZERO => Re { re: T::infinity() },
                            _ => One,
                        },
                        Im { im: i2 } => match i2 {
                            _ if im == T::ZERO => One,
                            _ => panic!("Zero to the complex power is undefined"),
                        },
                        ReIm { re, im: i2 } => match i2 {
                            _ if i2 == T::ZERO => match re {
                                _ if re > T::ZERO => Zero,
                                _ if re < T::ZERO => Re { re: T::infinity() },
                                _ => One,
                            },
                            _ => panic!("Zero to the complex power is undefined"),
                        },
                        Ph { ph } => match ph.principal_val() {
                            x if x == T::ZERO => Zero,
                            x if x == T::PI => Re { re: T::infinity() },
                            _ => panic!("Zero to the complex power is undefined"),
                        },
                        Pl { rad: rad2, ph } => match rad2 {
                            x if x == T::ZERO => One,
                            _ => match ph.principal_val() {
                                x if x == T::ZERO => Zero,
                                x if x == T::PI => Re { re: T::infinity() },
                                _ => panic!("Zero to the complex power is undefined"),
                            },
                        },
                        Ln { re, im: i2 } => match re {
                            x if x < T::LN_THRESHOLD => One,
                            _ => match i2.principal_val() {
                                x if x == T::ZERO => Zero,
                                x if x == T::PI => Re { re: T::infinity() },
                                _ => panic!("Zero to the complex power is undefined"),
                            },
                        },
                        _ => unreachable!(),
                    },
                    _ if rad == T::ONE => {
                        let ph = im.atan2(re);
                        match x {
                            J => Re { re: (-ph).exp() },
                            NegJ => Re { re: ph.exp() },
                            Re { re: r2 } => Ph { ph: ph * r2 },
                            Im { im: i2 } => Re {
                                re: (-ph * i2).exp(),
                            },
                            ReIm { re: r2, im: i2 } => Ln {
                                re: -ph * i2,
                                im: ph * r2,
                            },
                            Ph { ph } => {
                                let (r2, i2) = (ph.cos(), ph.sin());
                                Ln {
                                    re: -ph * i2,
                                    im: ph * r2,
                                }
                            }
                            Pl { rad: rad2, ph } => {
                                let (r2, i2) = (rad2 * ph.cos(), rad2 * ph.sin());
                                Ln {
                                    re: -ph * i2,
                                    im: ph * r2,
                                }
                            }
                            Ln {
                                re: ln_rad2,
                                im: ph,
                            } => {
                                let rad2 = ln_rad2.exp();
                                let (r2, i2) = (rad2 * ph.cos(), rad2 * ph.sin());
                                Ln {
                                    re: -ph * i2,
                                    im: ph * r2,
                                }
                            }
                            _ => unreachable!(),
                        }
                    }
                    _ => {
                        let ln_rad = rad.ln();
                        let ph = im.atan2(re);
                        match x {
                            J => Ln {
                                re: -ph,
                                im: ln_rad,
                            },
                            NegJ => Ln { re: ph, im: ln_rad },
                            Re { re: r2 } => Ln {
                                re: ln_rad * r2,
                                im: ph * r2,
                            },
                            Im { im: i2 } => Ln {
                                re: -ph * i2,
                                im: ln_rad * i2,
                            },
                            ReIm { re: r2, im: i2 } => Ln {
                                re: ln_rad * r2 - ph * i2,
                                im: ln_rad * i2 + ph * r2,
                            },
                            Ph { ph } => {
                                let (r2, i2) = (ph.cos(), ph.sin());
                                Ln {
                                    re: ln_rad * r2 - ph * i2,
                                    im: ln_rad * i2 + ph * r2,
                                }
                            }
                            Pl { rad: rad2, ph } => {
                                let (r2, i2) = (rad2 * ph.cos(), rad2 * ph.sin());
                                Ln {
                                    re: ln_rad * r2 - ph * i2,
                                    im: ln_rad * i2 + ph * r2,
                                }
                            }
                            Ln {
                                re: ln_rad2,
                                im: ph,
                            } => {
                                let rad2 = ln_rad2.exp();
                                let (r2, i2) = (rad2 * ph.cos(), rad2 * ph.sin());
                                Ln {
                                    re: ln_rad * r2 - ph * i2,
                                    im: ln_rad * i2 + ph * r2,
                                }
                            }
                            _ => unreachable!(),
                        }
                    }
                }
            }
        }
    }
}

macro_rules! impl_cpx_hash {
    ($float:ty) => {
        impl core::hash::Hash for Cpx<$float> {
            fn hash<H: core::hash::Hasher>(&self, state: &mut H) {
                use Cpx::*;
                let canonicalized = self.canonicalize_lazy();
                match canonicalized {
                    Zero => 0u8.hash(state),
                    One => 1u8.hash(state),
                    NegOne => 2u8.hash(state),
                    J => 3u8.hash(state),
                    NegJ => 4u8.hash(state),
                    Re { re } => {
                        5u8.hash(state);
                        re.to_bits().hash(state);
                    }
                    Im { im } => {
                        6u8.hash(state);
                        im.to_bits().hash(state);
                    }
                    ReIm { re, im } => {
                        7u8.hash(state);
                        re.to_bits().hash(state);
                        im.to_bits().hash(state);
                    }
                    Ph { ph } => {
                        8u8.hash(state);
                        ph.to_bits().hash(state);
                    }
                    Pl { rad, ph } => {
                        9u8.hash(state);
                        rad.to_bits().hash(state);
                        ph.to_bits().hash(state);
                    }
                    Ln { .. } => unreachable!(),
                }
            }
        }
    };
}

impl_cpx_hash!(f32);
impl_cpx_hash!(f64);

macro_rules! comm {
    (($a:pat, $b:pat)) => {
        ($a, $b) | ($b, $a)
    };
}
impl<T: FloatExt> PartialEq for Cpx<T> {
    /// # Examples
    ///
    /// ```
    /// # use cpx_coords::*;
    /// let a = Cpx::Zero;
    /// let b = Cpx::Zero;
    /// let c = Cpx::Re { re: 1.0 };
    /// let d = Cpx::Re { re: 1.0 };
    /// let e = Cpx::Im { im: 1.0 };
    /// let f = Cpx::Im { im: 2.0 };
    /// let g = Cpx::ReIm { re: 1.0, im: 0.0 };
    /// let h = Cpx::ReIm { re: 1.0, im: 0.0 };
    /// let i = Cpx::Pl { rad: 2.0, ph: 3.0 };
    /// let j = Cpx::Pl { rad: 2.0, ph: 3.0 };
    ///
    /// assert_eq!(a,b);
    /// assert_eq!(c,d);
    /// assert_eq!(g,h);
    /// assert_eq!(i,j);
    /// assert_ne!(e,f);
    /// assert_ne!(a,c);
    /// ```
    ///
    /// # Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // Identity cases:
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = zero == zero; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = neg_j == neg_j; }, || { let _ = neg_j; }, batch, reps)));
    /// // assert!((..=14).contains(&measure_cycles( || { let _ = re == re; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=8).contains(&measure_cycles( || { let _ = im == im; }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=18).contains(&measure_cycles( || { let _ = ph == ph; }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=25).contains(&measure_cycles( || { let _ = pl == pl; }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=25).contains(&measure_cycles( || { let _ = ln == ln; }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=16).contains(&measure_cycles( || { let _ = re_im == re_im; }, || { let _ = re_im; }, batch, reps)));
    ///
    /// // Crossing cases:
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = zero == one; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = zero == ph; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = zero == re; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=12).contains(&measure_cycles( || { let _ = zero == re_im; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = zero == pl; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = zero == ln; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = one == re; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=15).contains(&measure_cycles( || { let _ = one == ph; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = one == pl; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = re == im; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=12).contains(&measure_cycles( || { let _ = re == re_im; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=15).contains(&measure_cycles( || { let _ = re == ph; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=18).contains(&measure_cycles( || { let _ = re == pl; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=18).contains(&measure_cycles( || { let _ = re == ln; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = ph == pl; }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = ph == ln; }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=45).contains(&measure_cycles( || { let _ = pl == ln; }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=30).contains(&measure_cycles( || { let _ = re_im == ph; }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=30).contains(&measure_cycles( || { let _ = re_im == pl; }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=50).contains(&measure_cycles( || { let _ = re_im == ln; }, || { let _ = re_im; }, batch, reps)));
    /// ```
    #[inline(always)]
    fn eq(&self, other: &Self) -> bool {
        use Cpx::*;
        if matches!(
            (self, other),
            (Zero, Zero) | (One, One) | (NegOne, NegOne) | (J, J) | (NegJ, NegJ)
        ) {
            return true;
        }
        if matches!(
            (self, other),
            comm!((Zero, One))
                | comm!((Zero, NegOne))
                | comm!((Zero, J))
                | comm!((Zero, NegJ))
                | comm!((Zero, Ph { .. }))
                | comm!((One, NegOne))
                | comm!((One, J))
                | comm!((One, NegJ))
                | comm!((One, Im { .. }))
                | comm!((NegOne, J))
                | comm!((NegOne, NegJ))
                | comm!((NegOne, Im { .. }))
                | comm!((J, NegJ))
                | comm!((J, Re { .. }))
                | comm!((NegJ, Re { .. }))
        ) {
            return false;
        }

        match (*self, *other) {
            //(Zero, Zero) | (One, One) | (NegOne, NegOne) | (J, J) | (NegJ, NegJ) => true,
            (Re { re: v0 }, Re { re: v1 }) | (Im { im: v0 }, Im { im: v1 }) => v0 == v1,
            (ReIm { re: v0, im: v1 }, ReIm { re: v2, im: v3 }) => v0 == v2 && v1 == v3,
            (Ph { ph: v0 }, Ph { ph: v1 }) => (v0 - v1).principal_val() == T::ZERO,
            (Pl { rad: v0, ph: v1 }, Pl { rad: v2, ph: v3 }) => match (v0, v2) {
                _ if v0 == T::ZERO && v2 == T::ZERO => true,
                _ if v0 == v2 => (v1 - v3).principal_val() == T::ZERO,
                _ => false,
            },
            (Ln { re: v0, im: v1 }, Ln { re: v2, im: v3 }) => match (v0, v2) {
                _ if v0 < T::LN_THRESHOLD && v2 < T::LN_THRESHOLD => true,
                _ if v0 == v2 => (v1 - v3).principal_val() == T::ZERO,
                _ => false,
            },
            //comm!((Zero, One))
            //| comm!((Zero, NegOne))
            //| comm!((Zero, J))
            //| comm!((Zero, NegJ))
            //| comm!((Zero, Ph { .. }))
            //| comm!((One, NegOne))
            //| comm!((One, J))
            //| comm!((One, NegJ))
            //| comm!((One, Im { .. }))
            //| comm!((NegOne, J))
            //| comm!((NegOne, NegJ))
            //| comm!((NegOne, Im { .. }))
            //| comm!((J, NegJ))
            //| comm!((J, Re { .. }))
            //| comm!((NegJ, Re { .. })) => false,
            comm!((Zero, Re { re: val }))
            | comm!((Zero, Im { im: val }))
            | comm!((Zero, Pl { rad: val, .. })) => val == T::ZERO,
            comm!((Zero, Ln { re: val, .. })) => val < T::LN_THRESHOLD,
            comm!((Zero, ReIm { re, im })) => re == T::ZERO && im == T::ZERO,
            comm!((One, Re { re: val })) | comm!((J, Im { im: val })) => val == T::ONE,
            comm!((NegOne, Re { re: val })) | comm!((NegJ, Im { im: val })) => val == T::NEG_ONE,
            comm!((One, ReIm { re: v0, im: v1 })) | comm!((J, ReIm { re: v1, im: v0 })) => {
                v0 == T::ONE && v1 == T::ZERO
            }
            comm!((NegOne, ReIm { re: v0, im: v1 })) | comm!((NegJ, ReIm { re: v1, im: v0 })) => {
                v0 == T::NEG_ONE && v1 == T::ZERO
            }
            comm!((One, Ph { ph })) => ph.principal_val() == T::ZERO,
            comm!((NegOne, Ph { ph })) => ph.principal_val() == T::PI,
            comm!((J, Ph { ph })) => ph.principal_val() == T::FRAC_PI_2,
            comm!((NegJ, Ph { ph })) => ph.principal_val() == T::NEG_FRAC_PI_2,
            comm!((One, Pl { rad, ph })) => {
                if rad != T::ONE {
                    false
                } else {
                    ph.principal_val() == T::ZERO
                }
            }
            comm!((NegOne, Pl { rad, ph })) => {
                if rad != T::ONE {
                    false
                } else {
                    ph.principal_val() == T::PI
                }
            }
            comm!((J, Pl { rad, ph })) => {
                if rad != T::ONE {
                    false
                } else {
                    ph.principal_val() == T::FRAC_PI_2
                }
            }
            comm!((NegJ, Pl { rad, ph })) => {
                if rad != T::ONE {
                    false
                } else {
                    ph.principal_val() == T::NEG_FRAC_PI_2
                }
            }
            comm!((One, Ln { re, im })) => {
                if re != T::ZERO {
                    false
                } else {
                    im.principal_val() == T::ZERO
                }
            }
            comm!((NegOne, Ln { re, im })) => {
                if re != T::ZERO {
                    false
                } else {
                    im.principal_val() == T::PI
                }
            }
            comm!((J, Ln { re, im })) => {
                if re != T::ZERO {
                    false
                } else {
                    im.principal_val() == T::FRAC_PI_2
                }
            }
            comm!((NegJ, Ln { re, im })) => {
                if re != T::ZERO {
                    false
                } else {
                    im.principal_val() == T::NEG_FRAC_PI_2
                }
            }
            comm!((Re { re }, Im { im })) => re == T::ZERO && im == T::ZERO,
            comm!((Re { re: v0 }, ReIm { re: v1, im: v2 }))
            | comm!((Im { im: v0 }, ReIm { re: v2, im: v1 })) => {
                if v2 != T::ZERO {
                    false
                } else {
                    v0 == v1
                }
            }
            comm!((Re { re }, Ph { ph })) => match re {
                r if r == T::ONE => ph.principal_val() == T::ZERO,
                r if r == T::NEG_ONE => ph.principal_val() == T::PI,
                _ => false,
            },
            comm!((Im { im }, Ph { ph })) => match im {
                i if i == T::ONE => ph.principal_val() == T::FRAC_PI_2,
                i if i == T::NEG_ONE => ph.principal_val() == T::NEG_FRAC_PI_2,
                _ => false,
            },
            comm!((Re { re }, Pl { rad, ph })) => match ph.principal_val() {
                p if p == T::ZERO => re == rad,
                p if p == T::PI => re == -rad,
                _ => false,
            },
            comm!((Im { im }, Pl { rad, ph })) => match ph.principal_val() {
                p if p == T::FRAC_PI_2 => im == rad,
                p if p == T::NEG_FRAC_PI_2 => im == -rad,
                _ => false,
            },
            comm!((Re { re }, Ln { re: ln_rad, im: ph })) => match ph.principal_val() {
                p if p == T::ZERO => re == ln_rad.exp(),
                p if p == T::PI => re == -ln_rad.exp(),
                _ => false,
            },
            comm!((Im { im }, Ln { re: ln_rad, im: ph })) => match ph.principal_val() {
                p if p == T::FRAC_PI_2 => im == ln_rad.exp(),
                p if p == T::NEG_FRAC_PI_2 => im == -ln_rad.exp(),
                _ => false,
            },
            comm!((Ph { ph: v0 }, Pl { rad, ph: v1 })) => match rad {
                r if r == T::ONE => (v0 - v1).principal_val() == T::ZERO,
                _ => false,
            },
            comm!((Ph { ph: v0 }, Ln { re: ln_rad, im: v1 })) => match ln_rad {
                r if r == T::ZERO => (v0 - v1).principal_val() == T::ZERO,
                _ => false,
            },
            comm!((Pl { rad, ph: v0 }, Ln { re: ln_rad, im: v1 })) => {
                match (v0 - v1).principal_val() {
                    _ if rad == T::ZERO && ln_rad < T::LN_THRESHOLD => true,
                    p if p == T::ZERO => rad.ln() == ln_rad,
                    _ => false,
                }
            }
            comm!((ReIm { re, im }, Ph { ph })) => {
                if re.hypot(im) != T::ONE {
                    false
                } else {
                    im.atan2(re) == ph.principal_val()
                }
            }
            comm!((ReIm { re, im }, Pl { rad, ph })) => {
                if re.hypot(im) != rad {
                    false
                } else {
                    im.atan2(re) == ph.principal_val()
                }
            }
            comm!((ReIm { re, im }, Ln { re: ln_rad, im: ph })) => {
                if re.hypot(im) != ln_rad.exp() {
                    false
                } else {
                    im.atan2(re) == ph.principal_val()
                }
            }
            _ => unreachable!(), //comm!((ReIm { re, im }, x)) => {
                                 //    let (rhs_re, rhs_im) = x.re_im();
                                 //    re == rhs_re && im == rhs_im
                                 //}
        }
    }
}

impl<T: FloatExt> Neg for Cpx<T> {
    type Output = Self;
    /// # Examples
    ///
    /// ```
    /// use cpx_coords::*;
    /// let z: Cpx<f64> = Cpx::Zero;
    /// //assert_eq!(-z, Cpx::Zero);
    ///
    /// let o: Cpx<f64> = Cpx::One;
    /// //assert_eq!(-o, Cpx::NegOne);
    ///
    /// let n: Cpx<f64> = Cpx::NegOne;
    /// //assert_eq!(-n, Cpx::One);
    ///
    /// let j: Cpx<f64> = Cpx::J;
    /// //assert_eq!(-j, Cpx::NegJ);
    ///
    /// let neg_j: Cpx<f64> = Cpx::NegJ;
    /// //assert_eq!(-neg_j, Cpx::J);
    ///
    /// let r = Cpx::Re { re: 3.0f32 };
    /// //assert_eq!(-r, Cpx::Re { re: -3.0 });
    ///
    /// let i = Cpx::Im { im: 2.0f32 };
    /// //assert_eq!(-i, Cpx::Im { im: -2.0 });
    ///
    /// let p = Cpx::Ph { ph: f32::FRAC_PI_2 };
    /// assert_eq!(-p, Cpx::Ph { ph: -f32::FRAC_PI_2 });
    ///
    /// let c = Cpx::ReIm { re: 1.0f32, im: 2.0 };
    /// //assert_eq!(-c, Cpx::ReIm { re: -1.0, im: -2.0 });
    ///
    /// let l = Cpx::Ln { re: 1.0f32, im: 3.0 };
    /// //assert_eq!(-l, Cpx::Ln { re: 1.0, im: (3.0 + f32::PI).principal_val() });
    ///
    /// let pl = Cpx::Pl { rad: 2.0f32, ph: 0.5 };
    /// //assert_eq!(-pl, Cpx::Pl { rad: 2.0, ph: (0.5 + f32::PI).principal_val() });
    /// ```
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let (zero, one, neg_one, j, neg_j) = (C::Zero, C::One, C::NegOne, C::J, C::NegJ);
    ///
    /// let (val,val_2): (F,F) = (3.0,4.0);
    /// let (re, im, re_im) = (C::Re { re: val}, C::Im { im: val}, C::ReIm { re: val,  im: val_2});
    ///
    /// let phase: F = F::FRAC_PI_6;
    /// let (ph, pl, ln) = (C::Ph { ph: phase}, C::Pl { rad: val.abs(), ph: phase}, C::Ln { re: val, im: phase});
    ///
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = -zero; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=1).contains(&measure_cycles( || { let _ = -one; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = -j; }, || { let _ = j; }, batch, reps)));
    /// // assert!((..=7).contains(&measure_cycles( || { let _ = -re; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = -im; }, || { let _ = im; }, batch, reps)));
    /// // assert!((..=11).contains(&measure_cycles( || { let _ = -re_im; }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = -ph; }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = -pl; }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = -ln; }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    fn neg(self) -> Self {
        use Cpx::*;
        if let Zero = self {
            return Zero;
        }
        if let One = self {
            return NegOne;
        }
        if let NegOne = self {
            return One;
        }
        if let J = self {
            return NegJ;
        }
        if let NegJ = self {
            return J;
        }
        match self {
            Re { re } => Re { re: -re },
            Im { im } => Im { im: -im },
            ReIm { re, im } => ReIm { re: -re, im: -im },
            Ph { ph } => match ph {
                _ if ph > T::ZERO => Ph { ph: ph - T::PI },
                _ => Ph { ph: ph + T::PI },
            },
            Pl { rad, ph } => match ph {
                _ if ph > T::ZERO => Pl {
                    rad,
                    ph: ph - T::PI,
                },
                _ => Pl {
                    rad,
                    ph: ph + T::PI,
                },
            },
            Ln { re, im: ph } => match ph {
                _ if ph > T::ZERO => Ln { re, im: ph - T::PI },
                _ => Ln { re, im: ph + T::PI },
            },
            _ => unreachable!(),
        }
    }
}

impl<T> Add for Cpx<T>
where
    T: FloatExt,
{
    type Output = Self;

    /// # Examples
    ///
    /// ```
    /// # use cpx_coords::*;
    /// use Cpx::*;
    ///
    /// let z = Zero;
    /// let o = One;
    /// let n = NegOne;
    /// let j = J;
    /// let nj = NegJ;
    ///
    /// assert_eq!(z + o, o);
    /// assert_eq!(o + n, Zero);
    /// assert_eq!(j + nj, Zero);
    ///
    /// assert_eq!(o + j, ReIm { re: 1.0, im: 1.0 });
    /// assert_eq!(n + nj, ReIm { re: -1.0, im: -1.0 });
    ///
    /// let r1 = Re { re: 1.5 };
    /// let r2 = Re { re: 0.5 };
    /// assert_eq!(r1 + r2, Re { re: 2.0 });
    ///
    /// let i1 = Im { im: 2.0 };
    /// let i2 = Im { im: 3.0 };
    /// assert_eq!(i1 + i2, Im { im: 5.0 });
    ///
    /// let cc1 = ReIm { re: 1.0, im: 2.0 };
    /// let cc2 = ReIm { re: 3.0, im: -1.0 };
    /// assert_eq!(cc1 + cc2, ReIm { re: 4.0, im: 1.0 });
    /// ```
    ///
    /// ## Performance
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let zero = C::Zero;
    /// let one = C::One;
    /// let neg_one = C::NegOne;
    /// let j = C::J;
    /// let neg_j = C::NegJ;
    /// let re = C::Re { re: 2.0};
    /// let im = C::Im { im: 3.0};
    /// let ph = C::Ph { ph: F::FRAC_PI_4};
    /// let ph2 = C::Ph { ph: F::FRAC_PI_6};
    /// let re_im = C::ReIm { re: 3.0,  im: 4.0};
    /// let pl = C::Pl { rad: 2.0, ph: F::FRAC_PI_6};
    /// let ln = C::Ln { re: 1.0, im: F::FRAC_PI_6};
    ///
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = zero + zero; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = one + one; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = one + j; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=6).contains(&measure_cycles( || { let _ = one + im; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=12).contains(&measure_cycles( || { let _ = one + re; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=12).contains(&measure_cycles( || { let _ = one + re_im; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=12).contains(&measure_cycles( || { let _ = re + re; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=12).contains(&measure_cycles( || { let _ = re + re_im; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=15).contains(&measure_cycles( || { let _ = re_im + re_im; }, || { let _ = re_im; }, batch, reps)));
    /// // assert!((..=55).contains(&measure_cycles( || { let _ = ph + ph2; }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=40).contains(&measure_cycles( || { let _ = ph + one; }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=50).contains(&measure_cycles( || { let _ = ph + neg_one; }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=55).contains(&measure_cycles( || { let _ = ph + j; }, || { let _ = ph; }, batch, reps)));
    ///
    /// // assert!((..=90).contains(&measure_cycles( || { let _ = pl + one; }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=120).contains(&measure_cycles( || { let _ = ln + one; }, || { let _ = ln; }, batch, reps)));
    /// // assert!((..=180).contains(&measure_cycles( || { let _ = pl + pl; }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=280).contains(&measure_cycles( || { let _ = ln + ln; }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    fn add(self, other: Self) -> Self {
        use Cpx::*;

        match (self, other) {
            comm!((Zero, x)) => x,
            comm!((One, NegOne)) | comm!((J, NegJ)) => Zero,
            (One, One) => Re { re: T::TWO },
            comm!((One, J)) => ReIm {
                re: T::ONE,
                im: T::ONE,
            },
            comm!((One, NegJ)) => ReIm {
                re: T::ONE,
                im: T::NEG_ONE,
            },
            (NegOne, NegOne) => Re { re: -T::TWO },
            comm!((NegOne, J)) => ReIm {
                re: T::NEG_ONE,
                im: T::ONE,
            },
            comm!((NegOne, NegJ)) => ReIm {
                re: T::NEG_ONE,
                im: T::NEG_ONE,
            },
            (J, J) => Im { im: T::TWO },
            (NegJ, NegJ) => Im { im: -T::TWO },
            comm!((One, Re { re })) => Re { re: re + T::ONE },
            comm!((One, Im { im })) => ReIm { re: T::ONE, im },
            comm!((One, ReIm { re, im })) => ReIm {
                re: re + T::ONE,
                im,
            },
            comm!((NegOne, Re { re })) => Re {
                re: re + T::NEG_ONE,
            },
            comm!((NegOne, Im { im })) => ReIm { re: T::NEG_ONE, im },
            comm!((NegOne, ReIm { re, im })) => ReIm {
                re: re + T::NEG_ONE,
                im,
            },
            comm!((J, Re { re })) => ReIm { re, im: T::ONE },
            comm!((J, Im { im })) => Im { im: im + T::ONE },
            comm!((J, ReIm { re, im })) => ReIm {
                re,
                im: im + T::ONE,
            },
            comm!((NegJ, Re { re })) => ReIm { re, im: T::NEG_ONE },
            comm!((NegJ, Im { im })) => Im {
                im: im + T::NEG_ONE,
            },
            comm!((NegJ, ReIm { re, im })) => ReIm {
                re,
                im: im + T::NEG_ONE,
            },
            (Re { re }, Re { re: r2 }) => Re { re: re + r2 },
            (Im { im }, Im { im: i2 }) => Im { im: im + i2 },
            comm!((Re { re }, Im { im })) => ReIm { re, im },
            comm!((ReIm { re, im }, Re { re: r2 })) => ReIm { re: re + r2, im },
            comm!((ReIm { re, im }, Im { im: i2 })) => ReIm { re, im: im + i2 },
            (ReIm { re, im }, ReIm { re: r2, im: i2 }) => ReIm {
                re: re + r2,
                im: im + i2,
            },
            (Ph { ph: ph1 }, Ph { ph: ph2 }) => {
                // Optimized summation of two unit-phase values in polar form.
                // Instead of converting to Cartesian and back (~170 cycles),
                // this uses the identity:
                //   e^{iφ₁} + e^{iφ₂} = 2cos((φ₁ - φ₂)/2) · e^{i(φ₁ + φ₂)/2}
                // Reduces cost to ~55 cycles by requiring only one cosine call.
                let (diff, ph) = (
                    (ph1 - ph2).principal_val() / T::TWO,
                    (ph1 + ph2).principal_val() / T::TWO,
                );
                let rad = T::TWO * diff.cos();
                Pl { rad, ph }
            }
            comm!((Ph { ph: ph1 }, One)) => {
                let diff = ph1.principal_val() / T::TWO;
                let rad = T::TWO * diff.cos();
                Pl { rad, ph: diff } // In this case, average phase equals phase difference.
            }
            comm!((Ph { ph: ph1 }, NegOne)) => {
                let diff = (ph1 - T::PI).principal_val() / T::TWO;
                let rad = T::TWO * diff.cos();
                Pl { rad, ph: diff } // In this case, average phase equals phase difference.
            }
            comm!((Ph { ph: ph1 }, J)) => {
                let (diff, ph) = (
                    (ph1 + T::NEG_FRAC_PI_2).principal_val() / T::TWO,
                    (ph1 + T::FRAC_PI_2).principal_val() / T::TWO,
                );
                let rad = T::TWO * diff.cos();
                Pl { rad, ph }
            }
            comm!((Ph { ph: ph1 }, NegJ)) => {
                let (diff, ph) = (
                    (ph1 + T::FRAC_PI_2).principal_val() / T::TWO,
                    (ph1 + T::NEG_FRAC_PI_2).principal_val() / T::TWO,
                );
                let rad = T::TWO * diff.cos();
                Pl { rad, ph }
            }
            comm!((One, x)) => {
                let (re, im) = x.re_im();
                ReIm {
                    re: re + T::ONE,
                    im,
                }
            }
            comm!((NegOne, x)) => {
                let (re, im) = x.re_im();
                ReIm {
                    re: re + T::NEG_ONE,
                    im,
                }
            }
            comm!((J, x)) => {
                let (re, im) = x.re_im();
                ReIm {
                    re,
                    im: im + T::ONE,
                }
            }
            comm!((NegJ, x)) => {
                let (re, im) = x.re_im();
                ReIm {
                    re,
                    im: im + T::NEG_ONE,
                }
            }
            comm!((Re { re }, x)) => {
                let (r2, im) = x.re_im();
                ReIm { re: re + r2, im }
            }
            comm!((Im { im }, x)) => {
                let (re, i2) = x.re_im();
                ReIm { re, im: im + i2 }
            }
            comm!((ReIm { re, im }, x)) => {
                let (r2, i2) = x.re_im();
                ReIm {
                    re: re + r2,
                    im: im + i2,
                }
            }
            (x, y) => {
                let (re, im) = x.re_im();
                let (r2, i2) = y.re_im();
                ReIm {
                    re: re + r2,
                    im: im + i2,
                }
            }
        }
    }
}
impl<T> Sub for Cpx<T>
where
    T: FloatExt,
{
    type Output = Self;
    /// # Examples
    ///
    /// ```
    /// use cpx_coords::*;
    /// use Cpx::*;
    ///
    /// let a = Re { re: 2.0 };
    /// let b = Re { re: 1.5 };
    /// let c = Im { im: 3.0 };
    /// let d = Im { im: 0.5 };
    /// let z = Zero;
    ///
    /// assert_eq!(a - b, Re { re: 0.5 });
    /// assert_eq!(c - d, Im { im: 2.5 });
    ///
    /// let cc1 = ReIm { re: 3.0, im: 4.0 };
    /// let cc2 = ReIm { re: 1.0, im: 1.0 };
    /// assert_eq!(cc1 - cc2, ReIm { re: 2.0, im: 3.0 });
    ///
    /// assert_eq!(a - z, a);
    /// assert_eq!(z - a, Re { re: -2.0 });
    /// ```
    #[inline(always)]
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl<T> Mul for Cpx<T>
where
    T: FloatExt,
{
    type Output = Self;
    /// Implements multiplication for symbolic complex types [`Cpx<T>`].
    ///
    /// # Examples
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt};
    /// use Cpx::*;
    ///
    /// type F = f64;
    ///
    ///  let result = Zero * Re { re: 1.0};
    ///  assert_eq!(result, Zero);
    ///
    ///  let result = One * Im { im: -1.0};
    ///  assert_eq!(result, Im { im: -1.0});
    ///
    ///  let result = NegOne * Ph { ph: F::FRAC_PI_2 };
    ///  assert_eq!(result, -Ph { ph: F::FRAC_PI_2 });
    ///
    ///  let result = J::<F> * J;
    ///  assert_eq!(result, NegOne);
    ///
    ///  let result = J::<F> * NegJ;
    ///  assert_eq!(result, One);
    ///
    ///  let result = J * Re { re: 1.0};
    ///  assert_eq!(result, J);
    ///
    ///  let result = J * Re { re: -1.0};
    ///  assert_eq!(result, NegJ);
    ///
    ///  let result = J * Re { re: 3.0};
    ///  assert_eq!(result, Im { im: 3.0});
    ///
    ///  let result = J * Im { im: 1.0};
    ///  assert_eq!(result, NegOne);
    ///
    ///  let result = J * Im { im: -1.0};
    ///  assert_eq!(result, One);
    ///
    ///  let result = J * Im { im: 2.0};
    ///  assert_eq!(result, Re { re: -2.0});
    ///
    ///  let result = NegJ * Re { re: 1.0};
    ///  assert_eq!(result, NegJ);
    ///
    ///  let result = NegJ * Re { re: -1.0};
    ///  assert_eq!(result, J);
    ///
    ///  let result = NegJ * Re { re: 1.5};
    ///  assert_eq!(result, Im { im: -1.5});
    ///
    ///  let result = NegJ * Im { im: 1.0};
    ///  assert_eq!(result, One);
    ///
    ///  let result = NegJ * Im { im: -1.0};
    ///  assert_eq!(result, NegOne);
    ///
    ///  let result = NegJ * Im { im: 0.25};
    ///  assert_eq!(result, Re { re: 0.25});
    ///
    ///  let result = Re { re: 2.0} * Re { re: 3.0};
    ///  assert_eq!(result, Re { re: 6.0});
    ///
    ///  let result = Re { re: 1.0} * Im { im: 3.0};
    ///  assert_eq!(result, Im { im: 3.0});
    ///
    ///  let result = Im { im: 2.0} * Im { im: 4.0};
    ///  assert_eq!(result, Re { re: -8.0});
    ///
    ///  let result = Ph { ph: F::FRAC_PI_4} * Ph { ph: F::FRAC_PI_4};
    ///  assert_eq!(result, Ph { ph: F::FRAC_PI_2});
    ///
    ///  let result = Ph { ph: F::FRAC_PI_4} * J;
    ///  assert_eq!(result, Ph { ph: F::FRAC_PI_4 * 3.0});
    ///
    ///  let result = Ph { ph: F::FRAC_PI_4} * NegJ;
    ///  assert_eq!(result, Ph { ph: -F::FRAC_PI_4});
    ///
    ///  let result = Ph { ph: F::FRAC_PI_4} * Re { re: 2.0};
    ///  assert_eq!(result, Pl { rad: 2.0, ph: F::FRAC_PI_4});
    ///
    ///  let result = Ph { ph: F::FRAC_PI_4} * Im { im: -2.0};
    ///  assert_eq!(result, Pl { rad: 2.0, ph: -F::FRAC_PI_4});
    ///
    ///  let result = Ln { re: 2.0_f64.ln(), im: F::FRAC_PI_6} * Ln { re: 3.0_f64.ln(), im: F::FRAC_PI_6};
    ///  assert_eq!(result.rad(), 6.0_f64);
    ///  assert_eq!(result.ph(), F::FRAC_PI_3);
    ///
    ///  let result = Ln { re: 2.0_f64.ln(), im: F::FRAC_PI_4} * J;
    ///  assert_eq!(result.ph(), (F::FRAC_PI_4 + F::FRAC_PI_2).principal_val());
    ///
    ///  let result = Ln { re: 2.0_f64.ln(), im: F::FRAC_PI_4} * Re { re: -1.0};
    ///  assert_eq!(result.ph(), (F::FRAC_PI_4 + F::PI).principal_val());
    ///
    ///  let result = Pl { rad: 2.0, ph: F::FRAC_PI_4} * J;
    ///  assert_eq!(result.ph(), (F::FRAC_PI_4 + F::FRAC_PI_2).principal_val());
    ///
    ///  let result = Pl { rad: 3.0, ph: F::FRAC_PI_4} * Pl { rad: 4.0, ph: F::FRAC_PI_2};
    ///  assert_eq!(result.rad(), 12.0);
    ///  assert_eq!(result.ph(), (F::FRAC_PI_4 + F::FRAC_PI_2).principal_val());
    ///
    ///  let result = Pl { rad: 2.0, ph: F::FRAC_PI_4} * Ln { re: 2.0_f64.ln(), im: F::FRAC_PI_2};
    ///  assert_eq!(result.rad(), 4.0);
    ///  assert_eq!(result.ph(), (F::FRAC_PI_4 + F::FRAC_PI_2).principal_val());
    /// ```
    ///
    /// ## Performance Characteristics
    ///
    /// ### Per-Variant Multiplication Costs
    ///
    /// As summarized in the table below, the baseline cost of multiplying two general complex numbers in Cartesian form (`Ccs * Ccs`) is ≤40 CPU cycles.
    /// In contrast, multiplications involving the most common symbolic constants, `Zero`, `One`, `J`, `NegJ`, typically used in Pauli matrices,
    /// complete in ≤10 cycles due to constant folding and early-return optimizations.
    ///
    /// When extending to Hadamard-like operations involving `FRAC_1_SQRT_2`, the relevant cases (`Real`, `Imag`) remain efficient, with costs ≤20 cycles.
    /// This implies that multiplications involving normalized real or imaginary values are 2×–4× faster than the general Cartesian case.
    ///
    /// To support a broader class of normalized complex numbers, we include `Phase` and `PL` (polar form with non-unit modulus).
    /// Multiplications between `Phase`, `PL`, and the symbolic constants (`Zero`, `One`, `J`, `NegJ`, `Real`, `Imag`) also complete in ≤20 cycles,
    /// maintaining the 2× performance advantage over general Cartesian operations.
    ///
    /// The most expensive multiplication path, consuming ~70–75 cycles, arises from fallback cases such as `Phase * Ccs`, where computing the phase
    /// requires a call to `atan2`. Despite its cost, we retain the `Ccs` form for its essential role in supporting complex addition and numeric fallback.
    ///
    /// Additionally, we include the `Ln` coordinate form to represent complex logarithms. While `Ln * Ln` multiplication is comparable in cost
    /// to `PL * PL` (~20 cycles), its primary benefit lies in facilitating exponential and power computations. Specifically, raising a complex number
    /// to another (`Cpx::powc`) is implemented efficiently via the identity:
    ///
    /// ```text
    ///    base^exp = exp(ln(base) * exp)
    /// ```
    ///
    /// Thus, `Ln` serves as a useful intermediate that pairs naturally with `Ccs` in exponentiation logic.
    ///
    /// ### Benchmark Harness (Partial)
    ///
    /// ```rust
    /// use cpx_coords::{Cpx, FloatExt, measure_cycles};
    ///
    /// let batch = 100_000;
    /// let reps = 50;
    ///
    /// type F = f64;
    /// type C = Cpx<F>;
    ///
    /// let zero = C::Zero;
    /// let one = C::One;
    /// let neg_one = C::NegOne;
    /// let j = C::J;
    /// let neg_j = C::NegJ;
    /// let re = C::Re { re: 2.0};
    /// let im = C::Im { im: 3.0};
    /// let ph = C::Ph { ph: F::FRAC_PI_4};
    /// let ph2 = C::Ph { ph: F::FRAC_PI_6};
    /// let re_im = C::ReIm { re: 3.0,  im: 4.0};
    /// let pl = C::Pl { rad: 2.0, ph: F::FRAC_PI_6};
    /// let ln = C::Ln { re: 1.0, im: F::FRAC_PI_6};
    ///
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = zero * zero; }, || { let _ = zero; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = one * one; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=2).contains(&measure_cycles( || { let _ = one * j; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = j * one; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = one * im; }, || { let _ = one; }, batch, reps)));
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = one * re_im; }, || { let _ = one; }, batch, reps)));
    ///
    /// // assert!((..=3).contains(&measure_cycles( || { let _ = neg_one * one; }, || { let _ = neg_one; }, batch, reps)));
    ///  assert!((..=17).contains(&measure_cycles( || { let _ = neg_one * j; }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=5).contains(&measure_cycles( || { let _ = j * neg_one; }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=20).contains(&measure_cycles( || { let _ = neg_one * im; }, || { let _ = neg_one; }, batch, reps)));
    /// // assert!((..=28).contains(&measure_cycles( || { let _ = neg_one * re_im; }, || { let _ = neg_one; }, batch, reps)));
    ///
    /// // assert!((17..=20).contains(&measure_cycles( || { let _ = re * neg_one; }, || { let _ = -re; }, batch, reps)));
    /// // assert!((0..=5).contains(&measure_cycles( || { let _ = re * neg_one; }, || { let _ = -re; }, batch, reps)));
    /// // assert!((17..=20).contains(&measure_cycles( || { let _ = im * neg_one; }, || { let _ = -im; }, batch, reps)));
    /// // assert!((20..=25).contains(&measure_cycles( || { let _ = re_im * neg_one; }, || { let _ = -re_im; }, batch, reps)));
    /// // assert!((20..=25).contains(&measure_cycles( || { let _ = ph * neg_one; }, || { let _ = -ph; }, batch, reps)));
    /// // assert!((20..=25).contains(&measure_cycles( || { let _ = pl * neg_one; }, || { let _ = -pl; }, batch, reps)));
    ///
    /// // assert!((..=8).contains(&measure_cycles( || { let _ = j * j; }, || { let _ = j; }, batch, reps)));
    /// // assert!((..=10).contains(&measure_cycles( || { let _ = j * re; }, || { let _ = j; }, batch, reps)));
    ///
    /// // assert!((..=13).contains(&measure_cycles( || { let _ = re * re; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=18).contains(&measure_cycles( || { let _ = re * re_im; }, || { let _ = re; }, batch, reps)));
    /// // assert!((..=35).contains(&measure_cycles( || { let _ = re_im * re_im; }, || { let _ = re_im; }, batch, reps)));
    ///
    /// // assert!((..=15).contains(&measure_cycles( || { let _ = ph * ph2; }, || { let _ = ph; }, batch, reps)));
    /// // assert!((..=55).contains(&measure_cycles( || { let _ = ph * j; }, || { let _ = ph; }, batch, reps)));
    ///
    /// // assert!((..=180).contains(&measure_cycles( || { let _ = pl * pl; }, || { let _ = pl; }, batch, reps)));
    /// // assert!((..=230).contains(&measure_cycles( || { let _ = ln * ln; }, || { let _ = ln; }, batch, reps)));
    /// ```
    #[inline(always)]
    fn mul(self, other: Self) -> Self {
        use Cpx::*;
        // Fast path for common constants
        if matches!(self, Zero) || matches!(other, Zero) {
            return Zero;
        }
        if let One = self {
            return other;
        }
        if let One = other {
            return self;
        }
        if let NegOne = self {
            return -other;
        }
        if let NegOne = other {
            return -self;
        }

        match (self, other) {
            //comm!((Zero, _)) => Zero,
            //comm!((One, x)) => x,w
            //comm!((NegOne, Re { re })) => Re { re: -re},
            //comm!((NegOne, x)) => -x,
            comm!((J, NegJ)) => One,
            (J, J) | (NegJ, NegJ) => NegOne,
            comm!((J, Re { re })) => Im { im: re },
            comm!((J, Im { im })) => Re { re: -im },
            comm!((J, ReIm { re, im })) => ReIm { re: -im, im: re },
            comm!((NegJ, Re { re })) => Im { im: -re },
            comm!((NegJ, Im { im })) => Re { re: im },
            comm!((NegJ, ReIm { re, im })) => ReIm { re: im, im: -re },
            comm!((J, Ph { ph })) => Ph {
                ph: ph + T::FRAC_PI_2,
            },
            comm!((J, Pl { rad, ph })) => Pl {
                rad,
                ph: ph + T::FRAC_PI_2,
            },
            comm!((J, Ln { re, im })) => Ln {
                re,
                im: im + T::FRAC_PI_2,
            },
            comm!((NegJ, Ph { ph })) => Ph {
                ph: ph + T::NEG_FRAC_PI_2,
            },
            comm!((NegJ, Pl { rad, ph })) => Pl {
                rad,
                ph: ph + T::NEG_FRAC_PI_2,
            },
            comm!((NegJ, Ln { re, im })) => Ln {
                re,
                im: im + T::NEG_FRAC_PI_2,
            },
            (Re { re }, Re { re: r2 }) => Re { re: re * r2 },
            comm!((Re { re }, Im { im })) => Im { im: re * im },
            comm!((Re { re }, ReIm { re: r2, im })) => ReIm {
                re: re * r2,
                im: re * im,
            },
            comm!((Re { re }, Ph { ph })) => match re {
                _ if re < T::ZERO => Pl {
                    rad: -re,
                    ph: ph + T::PI,
                },
                _ => Pl { rad: re, ph },
            },
            comm!((Re { re }, Pl { rad, ph })) => match rad * re {
                x if x < T::ZERO => Pl {
                    rad: -x,
                    ph: ph + T::PI,
                },
                x => Pl { rad: x, ph },
            },
            comm!((Re { re: r2 }, Ln { re, im })) => match r2 * re.exp() {
                x if x < T::ZERO => Pl {
                    rad: -x,
                    ph: im + T::PI,
                },
                x => Pl { rad: x, ph: im },
            },
            (Im { im }, Im { im: i2 }) => Re { re: -im * i2 },
            comm!((Im { im }, ReIm { re, im: i2 })) => ReIm {
                re: -im * i2,
                im: im * re,
            },
            comm!((Im { im }, Ph { ph })) => match im {
                x if x < T::ZERO => Pl {
                    rad: -x,
                    ph: ph + T::NEG_FRAC_PI_2,
                },
                x => Pl {
                    rad: x,
                    ph: ph + T::FRAC_PI_2,
                },
            },
            comm!((Im { im }, Pl { rad, ph })) => match rad * im {
                x if x < T::ZERO => Pl {
                    rad: -x,
                    ph: ph + T::NEG_FRAC_PI_2,
                },
                x => Pl {
                    rad: x,
                    ph: ph + T::FRAC_PI_2,
                },
            },
            comm!((Im { im }, Ln { re, im: i2 })) => match im * re.exp() {
                x if x < T::ZERO => Pl {
                    rad: -x,
                    ph: i2 + T::NEG_FRAC_PI_2,
                },
                x => Pl {
                    rad: x,
                    ph: i2 + T::FRAC_PI_2,
                },
            },
            comm!((Ph { ph: ph2 }, Pl { rad, ph })) => Pl { rad, ph: ph + ph2 },
            comm!((Ph { ph }, Ln { re, im })) => Ln { re, im: im + ph },
            comm!((Ph { ph }, ReIm { re, im })) => Pl {
                rad: re.hypot(im),
                ph: ph + im.atan2(re),
            },
            (Ph { ph }, Ph { ph: p2 }) => Ph { ph: ph + p2 },
            (Pl { rad, ph }, Pl { rad: rad2, ph: ph2 }) => Pl {
                rad: rad * rad2,
                ph: ph + ph2,
            },
            (Ln { re, im }, Ln { re: r2, im: i2 }) => Ln {
                re: re + r2,
                im: im + i2,
            },
            comm!((Pl { rad, ph }, Ln { re, im })) => Pl {
                rad: rad * re.exp(),
                ph: ph + im,
            },
            comm!((Pl { rad, ph }, ReIm { re, im })) => Pl {
                rad: rad * re.hypot(im),
                ph: ph + im.atan2(re),
            },
            comm!((Ln { re, im }, ReIm { re: r2, im: i2 })) => Pl {
                rad: re.exp() + r2.hypot(i2),
                ph: im + i2.atan2(r2),
            },
            (ReIm { re, im }, ReIm { re: r2, im: i2 }) => ReIm {
                re: (re * r2) - (im * i2),
                im: (re * i2) + (im * r2),
            },
            _ => unreachable!(),
        }
    }
}

#[allow(clippy::suspicious_arithmetic_impl)]
impl<T> Div for Cpx<T>
where
    T: FloatExt,
{
    type Output = Self;

    fn div(self, other: Self) -> Self {
        self * other.recip()
    }
}

impl<T> Add<T> for Cpx<T>
where
    T: FloatExt,
{
    type Output = Self;
    fn add(self, other: T) -> Self {
        use Cpx::*;
        self + Re { re: other }
    }
}

impl<T> Sub<T> for Cpx<T>
where
    T: FloatExt,
{
    type Output = Self;
    fn sub(self, other: T) -> Self {
        self + (-other)
    }
}

impl<T> Mul<T> for Cpx<T>
where
    T: FloatExt,
{
    type Output = Self;
    fn mul(self, other: T) -> Self {
        self * Cpx::Re { re: other }
    }
}

impl<T> Div<T> for Cpx<T>
where
    T: FloatExt,
{
    type Output = Self;
    fn div(self, other: T) -> Self {
        self / Cpx::Re { re: other }
    }
}
macro_rules! impl_cpx_assign_ops {
    // Implements trait for Cpx<T> <op>= Cpx<T>
    ($trait:ident, $method:ident, $op:tt, Self) => {
        impl<T> core::ops::$trait for Cpx<T>
        where
            T: FloatExt,
        {
            fn $method(&mut self, rhs: Self) {
                *self = *self $op rhs;
            }
        }
    };

    // Implements trait for Cpx<T> <op>= T
    ($trait:ident, $method:ident, $op:tt, $rhs:ty) => {
        impl<T> core::ops::$trait<$rhs> for Cpx<T>
        where
            T: FloatExt,
        {
            fn $method(&mut self, rhs: $rhs) {
                *self = *self $op rhs;
            }
        }
    };
}

impl_cpx_assign_ops!(AddAssign, add_assign, +, Self);
impl_cpx_assign_ops!(SubAssign, sub_assign, -, Self);
impl_cpx_assign_ops!(MulAssign, mul_assign, *, Self);
impl_cpx_assign_ops!(DivAssign, div_assign, /, Self);

impl_cpx_assign_ops!(AddAssign, add_assign, +, T);
impl_cpx_assign_ops!(SubAssign, sub_assign, -, T);
impl_cpx_assign_ops!(MulAssign, mul_assign, *, T);
impl_cpx_assign_ops!(DivAssign, div_assign, /, T);

macro_rules! measure_batch {
    ($func:expr, $batch:expr) => {{
        unsafe { __cpuid(0) };
        let start = unsafe { _rdtsc() };
        for _ in 0..$batch {
            $func();
            black_box(());
        }
        let end = unsafe { __rdtscp(&mut 0) };
        unsafe { __cpuid(0) };
        end.saturating_sub(start)
    }};
}

macro_rules! measure_cycles_common {
    ($f:expr, $g:expr, $batch:expr, $reps:expr) => {{
        let mut total: u64 = 0;
        for _ in 0..$reps {
            let total_f = measure_batch!($f, $batch);
            let total_g = measure_batch!($g, $batch);
            let delta = total_f.saturating_sub(total_g);
            total += delta / $batch;
        }
        total / $reps as u64
    }};
}

#[inline(never)]
pub fn measure_cycles<F, G>(f: F, g: G, batch: u64, reps: u32) -> u64
where
    F: Fn(),
    G: Fn(),
{
    use core::arch::x86_64::{__cpuid, __rdtscp, _rdtsc};
    use core::hint::black_box;

    measure_cycles_common!(f, g, batch, reps)
}

#[inline(never)]
pub fn measure_cycles_mut<F, G>(mut f: F, mut g: G, batch: u64, reps: u32) -> u64
where
    F: FnMut(),
    G: FnMut(),
{
    use core::arch::x86_64::{__cpuid, __rdtscp, _rdtsc};
    use core::hint::black_box;

    measure_cycles_common!(f, g, batch, reps)
}

/// Private module to seal the `UnsignedIndex` trait.
mod sealed {
    pub trait Sealed {}
    impl Sealed for u32 {}
    impl Sealed for u64 {}
}

use core::convert::{TryFrom, TryInto};
use core::fmt::Debug;

pub trait CpxAsKey:
    sealed::Sealed
    + Copy
    + Ord
    + PartialOrd
    + Hash
    + Eq
    + PartialEq
    + TryFrom<usize>
    + TryInto<usize>
    + 'static
{
}

impl CpxAsKey for u32 {}
impl CpxAsKey for u64 {}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CpxKey<T: CpxAsKey> {
    Zero,
    One,
    NegOne,
    J,
    NegJ,
    Re { re_bits: T },
    Im { im_bits: T },
    Ph { ph_bits: T },
    ReIm { re_bits: T, im_bits: T },
    Ln { re_bits: T, im_bits: T },
    Pl { rad_bits: T, ph_bits: T },
}

macro_rules! impl_cpx_to_key {
    ($float:ty => $uint:ty) => {
        impl From<Cpx<$float>> for CpxKey<$uint> {
            fn from(z: Cpx<$float>) -> Self {
                match z {
                    Cpx::Zero => CpxKey::Zero,
                    Cpx::One => CpxKey::One,
                    Cpx::NegOne => CpxKey::NegOne,
                    Cpx::J => CpxKey::J,
                    Cpx::NegJ => CpxKey::NegJ,
                    Cpx::Re { re } => CpxKey::Re {
                        re_bits: re.to_bits(),
                    },
                    Cpx::Im { im } => CpxKey::Im {
                        im_bits: im.to_bits(),
                    },
                    Cpx::Ph { ph } => CpxKey::Ph {
                        ph_bits: ph.to_bits(),
                    },
                    Cpx::ReIm { re, im } => CpxKey::ReIm {
                        re_bits: re.to_bits(),
                        im_bits: im.to_bits(),
                    },
                    Cpx::Ln { re, im } => CpxKey::Ln {
                        re_bits: re.to_bits(),
                        im_bits: im.to_bits(),
                    },
                    Cpx::Pl { rad, ph } => CpxKey::Pl {
                        rad_bits: rad.to_bits(),
                        ph_bits: ph.to_bits(),
                    },
                }
            }
        }
    };
}

impl_cpx_to_key!(f32 => u32);
impl_cpx_to_key!(f64 => u64);

macro_rules! impl_key_to_cpx {
    ($uint:ty => $float:ty) => {
        impl From<CpxKey<$uint>> for Cpx<$float> {
            fn from(z: CpxKey<$uint>) -> Self {
                match z {
                    CpxKey::Zero => Cpx::Zero,
                    CpxKey::One => Cpx::One,
                    CpxKey::NegOne => Cpx::NegOne,
                    CpxKey::J => Cpx::J,
                    CpxKey::NegJ => Cpx::NegJ,
                    CpxKey::Re { re_bits } => Cpx::Re {
                        re: <$float>::from_bits(re_bits),
                    },
                    CpxKey::Im { im_bits } => Cpx::Im {
                        im: <$float>::from_bits(im_bits),
                    },
                    CpxKey::Ph { ph_bits } => Cpx::Ph {
                        ph: <$float>::from_bits(ph_bits),
                    },
                    CpxKey::ReIm { re_bits, im_bits } => Cpx::ReIm {
                        re: <$float>::from_bits(re_bits),
                        im: <$float>::from_bits(im_bits),
                    },
                    CpxKey::Ln { re_bits, im_bits } => Cpx::Ln {
                        re: <$float>::from_bits(re_bits),
                        im: <$float>::from_bits(im_bits),
                    },
                    CpxKey::Pl { rad_bits, ph_bits } => Cpx::Pl {
                        rad: <$float>::from_bits(rad_bits),
                        ph: <$float>::from_bits(ph_bits),
                    },
                }
            }
        }
    };
}

impl_key_to_cpx!(u32 => f32);
impl_key_to_cpx!(u64 => f64);

impl<T: CpxAsKey> CpxKey<T> {
    /// ````
    /// use cpx_coords::{Cpx, FloatExt, CpxKey, CpxAsKey};
    ///
    /// // A simple test that CpxKey can be treated as the Key of BTreeMap. Therefore, a map (CpxKey, separable_state) can be used to denote superposition state.
    /// /// use std::collections::BTreeMap;
    /// //let mut map: BTreeMap<CpxKey<u32>, i32> = BTreeMap::new();
    /// //map.insert(CpxKey::One, 42);
    ///
    /// // A simple test that CpxKey can be treated as the Key of HashMap. Therefore, a map (CpxKey, separable_state) can be used to denote superposition state.
    /// //use std::collections::HashMap;
    /// //let mut map: HashMap<CpxKey<u32>,i32> = HashMap::new();
    /// //map.insert(CpxKey::One, 42);
    ///
    /// type F = f64;
    /// let re_im: Cpx<F> = Cpx::ReIm { re: 0.6, im: 0.8};
    /// let key: CpxKey<u64> = CpxKey::from(re_im.clone());
    /// let recovered: Cpx<F> = Cpx::from(key);
    /// assert_eq!(re_im,recovered);
    /// ````
    pub fn null(self) -> Self {
        self
    }
}
