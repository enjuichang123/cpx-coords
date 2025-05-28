// Copyright 2025 En-Jui Chang
//
// Licensed under the Apache License, Version 2.0 <http://www.apache.org/licenses/LICENSE-2.0>
// or the MIT license <http://opensource.org/licenses/MIT>, at your option.
// This file may not be copied, modified, or distributed except according to those terms.

use core::f32::consts::{
    E as E32, FRAC_PI_2 as FRAC_PI_2_32, FRAC_PI_4 as FRAC_PI_4_32, PI as PI32,
    SQRT_2 as SQRT_2_32, TAU as TAU32,
};
use core::f64::consts::{
    E as E64, FRAC_PI_2 as FRAC_PI_2_64, FRAC_PI_4 as FRAC_PI_4_64, PI as PI64,
    SQRT_2 as SQRT_2_64, TAU as TAU64,
};
use core::hash::{Hash, Hasher};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};
use num_traits::{Float, NumCast};

/// Extension trait for providing common mathematical constants for floats.
pub trait FloatExt: Float + NumCast + Copy + 'static {
    /// Returns the value of π.
    fn pi() -> Self;
    /// Returns the value of τ (2π).
    fn tau() -> Self;
    /// Returns the base of natural logarithms, *e*.
    fn e() -> Self;
    /// Returns the square root of 2.
    fn sqrt2() -> Self;
    /// Returns π/2.
    fn frac_pi_2() -> Self;
    /// Returns π/4.
    fn frac_pi_4() -> Self;
    /// A helper function for branch cut.
    fn rem_euclid(a: Self, b: Self) -> Self {
        let r = a % b;
        if r < Self::zero() { r + b.abs() } else { r }
    }
    /// Wraps a phase angle into the (-π, π] range.
    fn branch_cut(ph: Self) -> Self {
        let wrapped = Self::rem_euclid(ph + Self::pi(), Self::tau()) - Self::pi();
        if wrapped <= -Self::pi() {
            wrapped + Self::tau()
        } else {
            wrapped
        }
    }
    /// Returns the threshold below which two floats are considered equal.
    fn threshold() -> Self;
    /// Approximately checks whether two `T` numbers are equal within a given threshold.
    ///
    /// # Arguments
    /// * `a` - First floating-point number.
    /// * `b` - Second floating-point number.
    ///
    /// # Returns
    /// * `true` if the absolute difference between `a` and `b` is less than or equal to `threshold`, otherwise `false`.
    fn approx_eq(a: Self, b: Self) -> bool {
        (a - b).abs() <= Self::threshold()
    }
}

impl FloatExt for f32 {
    fn pi() -> Self {
        PI32
    }
    fn tau() -> Self {
        TAU32
    }
    fn e() -> Self {
        E32
    }
    fn sqrt2() -> Self {
        SQRT_2_32
    }
    fn frac_pi_2() -> Self {
        FRAC_PI_2_32
    }
    fn frac_pi_4() -> Self {
        FRAC_PI_4_32
    }
    fn threshold() -> Self {
        1e-6f32
    }
}

impl FloatExt for f64 {
    fn pi() -> Self {
        PI64
    }
    fn tau() -> Self {
        TAU64
    }
    fn e() -> Self {
        E64
    }
    fn sqrt2() -> Self {
        SQRT_2_64
    }
    fn frac_pi_2() -> Self {
        FRAC_PI_2_64
    }
    fn frac_pi_4() -> Self {
        FRAC_PI_4_64
    }
    fn threshold() -> Self {
        1e-12f64
    }
}

/// Enum representing complex numbers in various coordinate systems.
#[derive(Debug, Clone, Copy)]
pub enum Cpx<T: FloatExt> {
    /// Represents the additive identity (0).
    Zero {},
    /// Represents the multiplicative identity (1).
    One {},
    /// Represents the negation of `One` (-1).
    NegOne {},
    /// Represents the imaginary unit (j), a square root of -1.
    J {},
    /// Represents the negation of `J` (-j).
    NegJ {},
    /// Represents a real number (re).
    Real { re: T },
    /// Represents a purely imaginary number (im * j).
    Imag { im: T },
    /// Represents a complex number with radius 1 and phase angle \(\phi \in (-\pi, \pi] \).
    Phase { ph: T },
    /// Represents a complex number in Cartesian coordinates (re + j * im).
    Ccs { re: T, im: T },
    /// Represents a complex number in logarithmic form: \(\ln(z) = re + j * im \),
    /// such that \(z = e^{re + j * im} \).
    Ln { re: T, im: T },
    /// Represents a complex number in polar coordinates: \(rad * e^{j * ph} \).
    PL { rad: T, ph: T },
}

/// Enum representing error types for complex number operations.
#[derive(Debug, PartialEq)]
pub enum CpxError {
    /// Occurs when attempting to divide by zero.
    DivisionByZero,
    // Additional error types can be added as needed.
}

impl<T: FloatExt> Hash for Cpx<T> {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let regularized = self.regularize();
        match regularized {
            Cpx::Zero {} => 0u8.hash(state),
            Cpx::One {} => 1u8.hash(state),
            Cpx::NegOne {} => 2u8.hash(state),
            Cpx::J {} => 3u8.hash(state),
            Cpx::NegJ {} => 4u8.hash(state),
            Cpx::Real { .. } => 5u8.hash(state),
            Cpx::Imag { .. } => 6u8.hash(state),
            Cpx::Phase { .. } => 7u8.hash(state),
            Cpx::Ccs { .. } => 8u8.hash(state),
            Cpx::Ln { .. } => 9u8.hash(state),
            Cpx::PL { .. } => 10u8.hash(state),
        }
    }
}

impl<T: FloatExt> PartialEq for Cpx<T> {
    fn eq(&self, other: &Self) -> bool {
        match (self.regularize(), other.regularize()) {
            (Cpx::Zero {}, Cpx::Zero {})
            | (Cpx::One {}, Cpx::One {})
            | (Cpx::NegOne {}, Cpx::NegOne {})
            | (Cpx::J {}, Cpx::J {})
            | (Cpx::NegJ {}, Cpx::NegJ {}) => true,
            (Cpx::Real { re: re1 }, Cpx::Real { re: re2 }) => re1 == re2,
            (Cpx::Imag { im: im1 }, Cpx::Imag { im: im2 }) => im1 == im2,
            (Cpx::Phase { ph: ph1 }, Cpx::Phase { ph: ph2 }) => ph1 == ph2,
            (Cpx::Ccs { re: re1, im: im1 }, Cpx::Ccs { re: re2, im: im2 })
            | (Cpx::Ln { re: re1, im: im1 }, Cpx::Ln { re: re2, im: im2 }) => {
                re1 == re2 && im1 == im2
            }
            (Cpx::PL { rad: rad1, ph: ph1 }, Cpx::PL { rad: rad2, ph: ph2 }) => {
                rad1 == rad2 && ph1 == ph2
            }
            _ => false,
        }
    }
}
impl<T: FloatExt> Eq for Cpx<T> {}

impl<T: FloatExt> From<(T, T)> for Cpx<T> {
    fn from((re, im): (T, T)) -> Self {
        Cpx::Ccs { re, im }.regularize()
    }
}

impl<T: FloatExt> Neg for Cpx<T> {
    type Output = Self;
    fn neg(self) -> Self {
        match self {
            Cpx::Zero {} => Cpx::Zero {},
            Cpx::One {} => Cpx::NegOne {},
            Cpx::NegOne {} => Cpx::One {},
            Cpx::J {} => Cpx::NegJ {},
            Cpx::NegJ {} => Cpx::J {},
            Cpx::Real { re } => Cpx::Real { re: -re },
            Cpx::Imag { im } => Cpx::Imag { im: -im },
            Cpx::Phase { ph } => Cpx::Phase {
                ph: T::branch_cut(ph + T::pi()),
            },
            Cpx::Ccs { re, im } => Cpx::Ccs { re: -re, im: -im },
            Cpx::Ln { re, im } => Cpx::Ln {
                re,
                im: T::branch_cut(im + T::pi()),
            },
            Cpx::PL { rad, ph } => Cpx::PL {
                rad,
                ph: T::branch_cut(ph + T::pi()),
            },
        }
    }
}

impl<T: FloatExt> Cpx<T> {
    /// Represents the complex number zero.
    pub const ZERO: Self = Cpx::Zero {};
    /// Represents the complex number one.
    pub const ONE: Self = Cpx::One {};
    /// Represents the imaginary unit \(j = \sqrt{-1} \).
    pub const J: Self = Cpx::J {};
    /// Represents the complex number negative one.
    pub const NEG_ONE: Self = Cpx::NegOne {};
    /// Represents the negative imaginary unit '-j'.
    pub const NEG_J: Self = Cpx::NegJ {};
    /// Represents the real number 1/sqrt(2).
    pub fn inv_sqrt_2() -> Self {
        Cpx::Real {
            re: T::one() / T::sqrt2(),
        }
    }
    /// Represents the square root of 'j'.
    pub fn sqrt_j() -> Self {
        Self::Phase { ph: T::frac_pi_4() }
    }

    /// Returns a canonicalized form of the complex number, reducing floating-point errors and mapping
    /// close values to standard representations (e.g., 0, 1, -1, j, -j).
    pub fn regularize(self) -> Self {
        match self {
            Cpx::Zero {} | Cpx::One {} | Cpx::NegOne {} | Cpx::J {} | Cpx::NegJ {} => self,
            Cpx::Real { re } => match () {
                _ if T::approx_eq(re, T::zero()) => Cpx::Zero {},
                _ if T::approx_eq(re, T::one()) => Cpx::One {},
                _ if T::approx_eq(re, -T::one()) => Cpx::NegOne {},
                _ => Cpx::Real { re },
            },
            Cpx::Imag { im } => match () {
                _ if T::approx_eq(im, T::zero()) => Cpx::Zero {},
                _ if T::approx_eq(im, T::one()) => Cpx::J {},
                _ if T::approx_eq(im, -T::one()) => Cpx::NegJ {},
                _ => Cpx::Imag { im },
            },
            Cpx::Phase { ph } => {
                let new_ph = T::branch_cut(ph);
                match () {
                    _ if T::approx_eq(new_ph, T::zero()) => Cpx::One {},
                    _ if T::approx_eq(new_ph, T::frac_pi_2()) => Cpx::J {},
                    _ if T::approx_eq(new_ph, T::pi()) => Cpx::NegOne {},
                    _ if T::approx_eq(new_ph, -T::frac_pi_2()) => Cpx::NegJ {},
                    _ => Cpx::Phase { ph: new_ph },
                }
            }
            Cpx::Ccs { re, im } => {
                let mag = re.hypot(im);
                let re_is_zero = T::approx_eq(re, T::zero());
                let im_is_zero = T::approx_eq(im, T::zero());

                match () {
                    _ if T::approx_eq(mag, T::zero()) => Cpx::Zero {},
                    _ if T::approx_eq(re, T::one()) && im_is_zero => Cpx::One {},
                    _ if T::approx_eq(re, -T::one()) && im_is_zero => Cpx::NegOne {},
                    _ if re_is_zero && T::approx_eq(im, T::one()) => Cpx::J {},
                    _ if re_is_zero && T::approx_eq(im, -T::one()) => Cpx::NegJ {},
                    _ if im_is_zero => Cpx::Real { re },
                    _ if re_is_zero => Cpx::Imag { im },
                    _ if T::approx_eq(mag, T::one()) => {
                        let phase = im.atan2(re);
                        Cpx::Phase { ph: phase }
                    }
                    _ => Cpx::Ccs { re, im },
                }
            }
            Cpx::Ln { re, im } => {
                let new_im = T::branch_cut(im);
                let re_is_zero = T::approx_eq(re, T::zero());

                match () {
                    _ if re < T::threshold().ln() => Cpx::Zero {},
                    _ if re_is_zero => match new_im {
                        x if T::approx_eq(x, T::zero()) => Cpx::One {},
                        x if T::approx_eq(x, T::frac_pi_2()) => Cpx::J {},
                        x if T::approx_eq(x, T::pi()) => Cpx::NegOne {},
                        x if T::approx_eq(x, -T::frac_pi_2()) => Cpx::NegJ {},
                        _ => Cpx::Phase { ph: new_im },
                    },
                    _ if T::approx_eq(new_im, T::zero()) => Cpx::Real { re: re.exp() },
                    _ if T::approx_eq(new_im, T::frac_pi_2()) => Cpx::Imag { im: re.exp() },
                    _ if T::approx_eq(new_im, T::pi()) => Cpx::Real { re: -re.exp() },
                    _ if T::approx_eq(new_im, -T::frac_pi_2()) => Cpx::Imag { im: -re.exp() },
                    _ => Cpx::Ln { re, im: new_im },
                }
            }
            Cpx::PL { rad, ph } => {
                let new_ph = T::branch_cut(ph);
                let rad_is_zero = T::approx_eq(rad, T::zero());
                let rad_is_one = T::approx_eq(rad, T::one());

                match () {
                    _ if rad_is_zero => Cpx::Zero {},
                    _ if rad_is_one && T::approx_eq(new_ph, T::zero()) => Cpx::One {},
                    _ if rad_is_one && T::approx_eq(new_ph, T::frac_pi_2()) => Cpx::J {},
                    _ if rad_is_one && T::approx_eq(new_ph, T::pi()) => Cpx::NegOne {},
                    _ if rad_is_one && T::approx_eq(new_ph, -T::frac_pi_2()) => Cpx::NegJ {},
                    _ if rad_is_one => Cpx::Phase { ph: new_ph },
                    _ if T::approx_eq(new_ph, T::zero()) => Cpx::Real { re: rad },
                    _ if T::approx_eq(new_ph, T::frac_pi_2()) => Cpx::Imag { im: rad },
                    _ if T::approx_eq(new_ph, T::pi()) => Cpx::Real { re: -rad },
                    _ if T::approx_eq(new_ph, -T::frac_pi_2()) => Cpx::Imag { im: -rad },
                    _ => Cpx::PL { rad, ph: new_ph },
                }
            }
        }
    }

    /// Returns the complex conjugate of the complex number.
    pub fn conj(self) -> Self {
        match self {
            Cpx::Zero {} | Cpx::One {} | Cpx::NegOne {} | Cpx::Real { .. } => self,
            Cpx::J {} => Cpx::NegJ {},
            Cpx::NegJ {} => Cpx::J {},
            Cpx::Imag { im } => Cpx::Imag { im: -im },
            Cpx::Phase { ph } => Cpx::Phase {
                ph: T::branch_cut(-ph),
            },
            Cpx::Ccs { re, im } => Cpx::Ccs { re, im: -im },
            Cpx::Ln { re, im } => Cpx::Ln {
                re,
                im: T::branch_cut(-im),
            },
            Cpx::PL { rad, ph } => Cpx::PL {
                rad,
                ph: T::branch_cut(-ph),
            },
        }
    }
    /// Returns the inverse (reciprocal) of the complex number.
    pub fn inv(&self) -> Result<Self, CpxError> {
        match self {
            Cpx::Zero {} => Err(CpxError::DivisionByZero),
            Cpx::One {} | Cpx::NegOne {} => Ok(*self),
            Cpx::J {} => Ok(Cpx::NegJ {}),
            Cpx::NegJ {} => Ok(Cpx::J {}),
            Cpx::Ccs { .. } => Ok(Cpx::PL {
                rad: T::one() / self.rad(),
                ph: T::branch_cut(-self.ph()),
            }),
            Cpx::Real { re } => Ok(Cpx::Real { re: T::one() / *re }),
            Cpx::Imag { im } => Ok(Cpx::Imag {
                im: -T::one() / *im,
            }),
            Cpx::Phase { ph } => Ok(Cpx::Phase {
                ph: T::branch_cut(-*ph),
            }),
            Cpx::Ln { re, im } => Ok(Cpx::Ln {
                re: -*re,
                im: T::branch_cut(-*im),
            }),
            Cpx::PL { rad, ph } => Ok(Cpx::PL {
                rad: T::one() / *rad,
                ph: T::branch_cut(-*ph),
            }),
        }
    }
    /// Returns the exponential of the complex number.
    pub fn exp(&self) -> Self {
        match self {
            Cpx::Zero {} => Cpx::One {},
            Cpx::One {} => Cpx::Real { re: T::e() },
            Cpx::NegOne {} => Cpx::Real {
                re: T::one() / T::e(),
            },
            Cpx::J {} => Cpx::Phase { ph: T::one() },
            Cpx::NegJ {} => Cpx::Phase { ph: -T::one() },
            Cpx::Real { re } => Cpx::Real { re: re.exp() },
            Cpx::Imag { im } => Cpx::Phase {
                ph: T::branch_cut(*im),
            }
            .regularize(),
            Cpx::Phase { ph } => Cpx::Ln {
                re: ph.cos(),
                im: ph.sin(),
            }
            .regularize(),
            Cpx::Ccs { re, im } => Cpx::Ln {
                re: *re,
                im: T::branch_cut(*im),
            },
            Cpx::Ln { re, im } => Cpx::Ln {
                re: re.exp() * im.cos(),
                im: T::branch_cut(re.exp() * im.sin()),
            }
            .regularize(),
            Cpx::PL { rad, ph } => Cpx::Ln {
                re: *rad * ph.cos(),
                im: T::branch_cut(*rad * ph.sin()),
            }
            .regularize(),
        }
    }
    /// Returns the logarithm of the complex number.
    pub fn ln(&self) -> Self {
        match self {
            Cpx::Zero {} => Cpx::Real {
                re: T::neg_infinity(),
            },
            Cpx::One {} => Cpx::Zero {},
            Cpx::NegOne {} => Cpx::Imag { im: T::pi() },
            Cpx::J {} => Cpx::Imag { im: T::frac_pi_2() },
            Cpx::NegJ {} => Cpx::Imag {
                im: -T::frac_pi_2(),
            },
            Cpx::Real { re } => {
                if *re >= T::zero() {
                    Cpx::Real { re: re.ln() }
                } else {
                    Cpx::Ccs {
                        re: (-*re).ln(),
                        im: T::pi(),
                    }
                }
            }
            Cpx::Imag { im } => {
                if *im >= T::zero() {
                    Cpx::Ccs {
                        re: im.ln(),
                        im: T::frac_pi_2(),
                    }
                } else {
                    Cpx::Ccs {
                        re: (-*im).ln(),
                        im: -T::frac_pi_2(),
                    }
                }
            }
            Cpx::Phase { ph } => Cpx::Imag { im: *ph },
            Cpx::Ccs { .. } => Cpx::Ccs {
                re: self.rad().ln(),
                im: self.ph(),
            },
            Cpx::Ln { re, im } => Cpx::Ccs { re: *re, im: *im },
            Cpx::PL { rad, ph } => Cpx::Ccs {
                re: rad.ln(),
                im: *ph,
            },
        }
    }
    /// Raises the complex number to an integer power using polar form.
    ///
    /// Computes `r^n * e^(i * n * θ)` by raising the magnitude to `n` and scaling the phase.
    /// Returns the regularized result.
    pub fn powi(self, n: i32) -> Self {
        use Cpx::*;

        if n == 0 {
            return One {};
        }

        let rad = self.rad().powi(n);
        let ph = self.ph() * T::from(n).expect("i32 to T conversion failed");
        PL { rad, ph }.regularize()
    }
    /// Raises the complex number to a floating-point power using polar form.
    ///
    /// Computes `r^f * e^(i * f * θ)` by raising the magnitude to `f` and scaling the phase.
    /// Returns the regularized result.
    pub fn powf(self, f: T) -> Self {
        use Cpx::*;

        if T::approx_eq(f, T::zero()) {
            return One {};
        }

        let rad = self.rad().powf(f);
        let ph = self.ph() * f;
        PL { rad, ph }.regularize()
    }
    /// Raises the complex number to a Cpx power using polar form.
    ///
    /// Computes `r^c * e^(i * c * θ)` by raising the magnitude to `c` and scaling the phase.
    /// Returns the regularized result.
    pub fn powc(self, c: Self) -> Self {
        use Cpx::*;

        // Anything raised to 0 is 1
        if matches!(c, Zero {}) {
            return One {};
        }

        // Zero raised to any nonzero is still zero (except negative powers, which you may want to check)
        if matches!(self, Zero {}) {
            return Zero {};
        }

        // Convert to polar logarithmic form: ln(r) + iθ
        let ln_self = Ln {
            re: self.rad().ln(),
            im: self.ph(),
        };

        // Exponentiation: e^{c * ln(self)}
        let ln_pow = ln_self * c;

        // Return exp(c * ln(self)), regularized
        Ln {
            re: ln_pow.re(),
            im: ln_pow.im(),
        }
        .regularize()
    }
    /// Returns the n-th root of the complex number.
    pub fn root_i(&self, n: i32) -> Self {
        self.powf(T::one() / T::from(n).unwrap())
    }
    /// Returns the real part of the complex number as T.
    pub fn re(&self) -> T {
        match self {
            Cpx::Zero {} | Cpx::J {} | Cpx::NegJ {} | Cpx::Imag { .. } => T::zero(),
            Cpx::One {} => T::one(),
            Cpx::NegOne {} => -T::one(),
            Cpx::Real { re } | Cpx::Ccs { re, .. } => *re,
            Cpx::Phase { ph } => ph.cos(),
            Cpx::Ln { re, im } => re.exp() * im.cos(),
            Cpx::PL { rad, ph } => *rad * ph.cos(),
        }
    }

    /// Returns the imaginary part of the complex number as `T`.
    pub fn im(&self) -> T {
        match self {
            Cpx::Zero {} | Cpx::One {} | Cpx::NegOne {} | Cpx::Real { .. } => T::zero(),
            Cpx::J {} => T::one(),
            Cpx::NegJ {} => -T::one(),
            Cpx::Imag { im } | Cpx::Ccs { im, .. } => *im,
            Cpx::Phase { ph } => ph.sin(),
            Cpx::Ln { re, im } => re.exp() * im.sin(),
            Cpx::PL { rad, ph } => *rad * ph.sin(),
        }
    }

    /// Returns the radius (magnitude) of the complex number as `T`.
    pub fn rad(&self) -> T {
        match self {
            Cpx::Zero {} => T::zero(),
            Cpx::One {} | Cpx::NegOne {} | Cpx::J {} | Cpx::NegJ {} | Cpx::Phase { .. } => T::one(),
            Cpx::Real { re } => re.abs(),
            Cpx::Imag { im } => im.abs(),
            Cpx::Ccs { re, im } => re.hypot(*im),
            Cpx::Ln { re, .. } => re.exp(),
            Cpx::PL { rad, .. } => *rad,
        }
    }

    /// Returns the phase angle of the complex number in the interval (-π, π] as `T`.
    pub fn ph(&self) -> T {
        match self {
            Cpx::Zero {} | Cpx::One {} => T::zero(),
            Cpx::J {} => T::frac_pi_2(),
            Cpx::NegOne {} => T::pi(),
            Cpx::NegJ {} => -T::frac_pi_2(),
            Cpx::Real { re } => {
                if *re >= T::zero() {
                    T::zero()
                } else {
                    T::pi()
                }
            }
            Cpx::Imag { im } => {
                if *im >= T::zero() {
                    T::frac_pi_2()
                } else {
                    -T::frac_pi_2()
                }
            }
            Cpx::Phase { ph } | Cpx::PL { ph, .. } => T::branch_cut(*ph),
            Cpx::Ln { im, .. } => T::branch_cut(*im),
            Cpx::Ccs { re, im } => im.atan2(*re),
        }
    }
    /// Factors the complex number into a radius and phase part.
    pub fn factor_out(&self) -> (Self, Self) {
        let re = self.rad();
        let ph = self.ph();
        (Self::Real { re }, Self::Phase { ph })
    }
}

impl<T> Add for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    type Output = Self;

    fn add(self, other: Self) -> Self {
        use Cpx::*;

        match (self, other) {
            (Zero {}, x) | (x, Zero {}) => x,
            (One {}, One {}) => Real {
                re: T::from(2).unwrap(),
            },
            (One {}, NegOne {}) | (NegOne {}, One {}) => Zero {},
            (One {}, J {}) | (J {}, One {}) => Ccs {
                re: T::one(),
                im: T::one(),
            },
            (One {}, NegJ {}) | (NegJ {}, One {}) => Ccs {
                re: T::one(),
                im: -T::one(),
            },
            (NegOne {}, NegOne {}) => Real {
                re: -T::from(2).unwrap(),
            },
            (NegOne {}, J {}) | (J {}, NegOne {}) => Ccs {
                re: -T::one(),
                im: T::one(),
            },
            (NegOne {}, NegJ {}) | (NegJ {}, NegOne {}) => Ccs {
                re: -T::one(),
                im: -T::one(),
            },
            (J {}, J {}) => Imag {
                im: T::from(2).unwrap(),
            },
            (J {}, NegJ {}) | (NegJ {}, J {}) => Zero {},
            (NegJ {}, NegJ {}) => Imag {
                im: -T::from(2).unwrap(),
            },
            (Real { re: r1 }, Real { re: r2 }) => Real { re: r1 + r2 },
            (Real { re }, One {}) | (One {}, Real { re }) => Real { re: re + T::one() },
            (Real { re }, NegOne {}) | (NegOne {}, Real { re }) => Real { re: re - T::one() },
            (Real { re }, J {}) | (J {}, Real { re }) => Ccs { re, im: T::one() },
            (Real { re }, NegJ {}) | (NegJ {}, Real { re }) => Ccs { re, im: -T::one() },
            (Imag { im: i1 }, Imag { im: i2 }) => Imag { im: i1 + i2 },
            (Real { re }, Imag { im }) | (Imag { im }, Real { re }) => Ccs { re, im },
            (One {}, Imag { im }) | (Imag { im }, One {}) => Ccs { re: T::one(), im },
            (NegOne {}, Imag { im }) | (Imag { im }, NegOne {}) => Ccs { re: -T::one(), im },
            (J {}, Imag { im }) | (Imag { im }, J {}) => Imag { im: im + T::one() },
            (NegJ {}, Imag { im }) | (Imag { im }, NegJ {}) => Imag { im: im - T::one() },
            (Ccs { re: r1, im: i1 }, Ccs { re: r2, im: i2 }) => Ccs {
                re: r1 + r2,
                im: i1 + i2,
            },
            (Ccs { re, im }, One {}) | (One {}, Ccs { re, im }) => Ccs {
                re: re + T::one(),
                im,
            },
            (Ccs { re, im }, NegOne {}) | (NegOne {}, Ccs { re, im }) => Ccs {
                re: re - T::one(),
                im,
            },
            (Ccs { re, im }, J {}) | (J {}, Ccs { re, im }) => Ccs {
                re,
                im: im + T::one(),
            },
            (Ccs { re, im }, NegJ {}) | (NegJ {}, Ccs { re, im }) => Ccs {
                re,
                im: im - T::one(),
            },
            (Ccs { re, im }, Real { re: r2 }) | (Real { re: r2 }, Ccs { re, im }) => {
                Ccs { re: re + r2, im }
            }
            (Ccs { re, im }, Imag { im: i2 }) | (Imag { im: i2 }, Ccs { re, im }) => {
                Ccs { re, im: im + i2 }
            }
            _ => Ccs {
                re: self.re() + other.re(),
                im: self.im() + other.im(),
            },
        }
    }
}
impl<T> AddAssign for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    fn add_assign(&mut self, rhs: Self) {
        *self = (*self) + rhs;
    }
}
impl<T> Sub for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

impl<T> SubAssign for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}

impl<T> Mul for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        match (self, other) {
            (Cpx::Zero {}, _) | (_, Cpx::Zero {}) => Cpx::Zero {}, // Anything * 0 = 0
            (Cpx::One {}, x) | (x, Cpx::One {}) => x,              // 1 * x = x
            (Cpx::NegOne {}, x) | (x, Cpx::NegOne {}) => -x,       // -1 * x = -x
            (Cpx::J {}, Cpx::J {}) => Cpx::NegOne {},              // j * j = -1
            (Cpx::J {}, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::J {}) => Cpx::One {}, // j * -j = 1
            (Cpx::NegJ {}, Cpx::NegJ {}) => Cpx::NegOne {},        // (-j) * (-j) = -1
            (Cpx::J {}, Cpx::Real { re }) | (Cpx::Real { re }, Cpx::J {}) => Cpx::Imag { im: re },
            (Cpx::J {}, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::J {}) => Cpx::Real { re: -im },
            (Cpx::NegJ {}, Cpx::Real { re }) | (Cpx::Real { re }, Cpx::NegJ {}) => {
                Cpx::Imag { im: -re }
            }
            (Cpx::NegJ {}, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::NegJ {}) => {
                Cpx::Real { re: im }
            }
            (Cpx::Real { re: r1 }, Cpx::Real { re: r2 }) => Cpx::Real { re: r1 * r2 }, // Real * Real
            (Cpx::Real { re }, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::Real { re }) => {
                Cpx::Imag { im: re * im } // Real * Imag or Imag * Real
            }
            (Cpx::Imag { im: i1 }, Cpx::Imag { im: i2 }) => Cpx::Real { re: -i1 * i2 }, // Imag * Imag = -Real
            (Cpx::Phase { ph }, Cpx::Phase { ph: p2 }) => Cpx::Phase { ph: ph + p2 },
            (Cpx::Phase { ph }, Cpx::J {}) | (Cpx::J {}, Cpx::Phase { ph }) => {
                Cpx::Phase { ph: ph + T::pi() }.regularize()
            }
            (Cpx::Phase { ph }, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::Phase { ph }) => {
                Cpx::Phase { ph: ph - T::pi() }.regularize()
            }
            (Cpx::Phase { ph }, Cpx::Real { re }) | (Cpx::Real { re }, Cpx::Phase { ph }) => {
                Cpx::PL { rad: re, ph }.regularize()
            }
            (Cpx::Phase { ph }, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::Phase { ph }) => {
                if im >= T::zero() {
                    Cpx::PL {
                        rad: im.abs(),
                        ph: ph + T::frac_pi_2(),
                    }
                    .regularize()
                } else {
                    Cpx::PL {
                        rad: im.abs(),
                        ph: ph - T::frac_pi_2(),
                    }
                    .regularize()
                }
            }
            (Cpx::Ln { re, im }, Cpx::Ln { re: r2, im: i2 }) => Cpx::Ln {
                re: re + r2,
                im: im + i2,
            },
            (Cpx::Ln { re, im }, Cpx::J {}) | (Cpx::J {}, Cpx::Ln { re, im }) => Cpx::Ln {
                re,
                im: im + T::frac_pi_2(),
            }
            .regularize(),
            (Cpx::Ln { re, im }, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::Ln { re, im }) => Cpx::Ln {
                re,
                im: im - T::frac_pi_2(),
            }
            .regularize(),
            (Cpx::Ln { re, im }, Cpx::Real { re: r2 })
            | (Cpx::Real { re: r2 }, Cpx::Ln { re, im }) => Cpx::Ln {
                re: re + r2.ln(),
                im,
            }
            .regularize(),
            (Cpx::Ln { re, im }, Cpx::Imag { im: i2 })
            | (Cpx::Imag { im: i2 }, Cpx::Ln { re, im }) => {
                if im >= T::zero() {
                    Cpx::Ln {
                        re: re + i2.abs().ln(),
                        im: im + T::frac_pi_2(),
                    }
                    .regularize()
                } else {
                    Cpx::Ln {
                        re: re + i2.abs().ln(),
                        im: im - T::frac_pi_2(),
                    }
                    .regularize()
                }
            }
            (Cpx::Ln { re, im }, Cpx::Phase { ph }) | (Cpx::Phase { ph }, Cpx::Ln { re, im }) => {
                Cpx::Ln { re, im: im + ph }.regularize()
            }
            (Cpx::PL { rad, ph }, Cpx::PL { rad: rad2, ph: ph2 }) => Cpx::PL {
                rad: rad * rad2,
                ph: ph + ph2,
            }
            .regularize(),
            (Cpx::PL { rad, ph }, Cpx::J {}) | (Cpx::J {}, Cpx::PL { rad, ph }) => Cpx::PL {
                rad,
                ph: ph + T::frac_pi_2(),
            }
            .regularize(),
            (Cpx::PL { rad, ph }, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::PL { rad, ph }) => Cpx::PL {
                rad,
                ph: ph - T::frac_pi_2(),
            }
            .regularize(),
            (Cpx::PL { rad, ph }, Cpx::Real { re }) | (Cpx::Real { re }, Cpx::PL { rad, ph }) => {
                Cpx::PL { rad: rad * re, ph }.regularize()
            }
            (Cpx::PL { rad, ph }, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::PL { rad, ph }) => {
                if im >= T::zero() {
                    Cpx::PL {
                        rad: rad * im.abs(),
                        ph: ph + T::frac_pi_2(),
                    }
                    .regularize()
                } else {
                    Cpx::PL {
                        rad: rad * im.abs(),
                        ph: ph - T::frac_pi_2(),
                    }
                    .regularize()
                }
            }
            (Cpx::PL { rad, ph }, Cpx::Phase { ph: ph2 })
            | (Cpx::Phase { ph: ph2 }, Cpx::PL { rad, ph }) => {
                Cpx::PL { rad, ph: ph + ph2 }.regularize()
            }
            (Cpx::PL { rad, ph }, Cpx::Ln { re, im })
            | (Cpx::Ln { re, im }, Cpx::PL { rad, ph }) => Cpx::Ln {
                re: re + rad.ln(),
                im: im + ph,
            }
            .regularize(),
            _ => Cpx::PL {
                rad: self.rad() * other.rad(),
                ph: self.ph() + other.ph(),
            }
            .regularize(),
        }
    }
}
impl<T> MulAssign for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}
impl<T> Div for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let new_rad = self.rad() / other.rad();
        let new_phase = T::branch_cut(self.ph() - other.ph());
        Cpx::PL {
            rad: new_rad,
            ph: new_phase,
        }
        .regularize()
    }
}

impl<T> DivAssign for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}

impl<T> Add<T> for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    type Output = Self;
    fn add(self, other: T) -> Self {
        self + Cpx::Real { re: other }
    }
}

impl<T> AddAssign<T> for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    fn add_assign(&mut self, rhs: T) {
        *self = *self + rhs;
    }
}
impl<T> Sub<T> for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    type Output = Self;
    fn sub(self, other: T) -> Self {
        self + (-other)
    }
}
impl<T> SubAssign<T> for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    fn sub_assign(&mut self, rhs: T) {
        *self = *self - rhs;
    }
}
impl<T> Mul<T> for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    type Output = Self;
    fn mul(self, other: T) -> Self {
        self * Cpx::Real { re: other }
    }
}
impl<T> MulAssign<T> for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    fn mul_assign(&mut self, rhs: T) {
        *self = *self * rhs;
    }
}
impl<T> Div<T> for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    type Output = Self;
    fn div(self, other: T) -> Self {
        self * (T::one() / other)
    }
}
impl<T> DivAssign<T> for Cpx<T>
where
    T: Float + NumCast + Copy + FloatExt,
{
    fn div_assign(&mut self, rhs: T) {
        *self = *self / rhs;
    }
}
