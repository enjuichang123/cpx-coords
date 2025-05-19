// Copyright 2025 En-Jui Chang
//
// Licensed under the Apache License, Version 2.0 <http://www.apache.org/licenses/LICENSE-2.0>
// or the MIT license <http://opensource.org/licenses/MIT>, at your option.
// This file may not be copied, modified, or distributed except according to those terms.

use core::f32::consts::{E, FRAC_PI_2, FRAC_PI_4, PI, SQRT_2, TAU};
use core::hash::{Hash, Hasher};
use core::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

/// Represents the complex number zero.
pub const ZERO: Cpx = Cpx::Zero {};
/// Represents the complex number one.
pub const ONE: Cpx = Cpx::One {};
/// Represents the complex number negative one.
pub const NEG_ONE: Cpx = Cpx::NegOne {};
/// Represents the imaginary unit \(j = \sqrt{-1} \).
pub const J: Cpx = Cpx::J {};
/// Represents the negative imaginary unit '-j'.
pub const NEG_J: Cpx = Cpx::NegJ {};
/// Represents the real number 1/sqrt(2).
pub const INV_SQRT_2: Cpx = Cpx::Real { re: 1.0 / SQRT_2 };
/// Represents the real number -1/sqrt(2).
pub const NEG_INV_SQRT_2: Cpx = Cpx::Real { re: -1.0 / SQRT_2 };
/// Represents the imaginary number j/sqrt(2).
pub const J_INV_SQRT_2: Cpx = Cpx::Imag { im: 1.0 / SQRT_2 };
/// Represents the imaginary number -j/sqrt(2).
pub const NEG_J_INV_SQRT_2: Cpx = Cpx::Imag { im: -1.0 / SQRT_2 };
/// Represents the square root of 'j'.
pub const SQRT_J: Cpx = Cpx::Phase { ph: FRAC_PI_4 };

/// Enum representing complex numbers in various coordinate systems.
#[derive(Debug, Clone, Copy)]
pub enum Cpx {
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
    Real { re: f32 },
    /// Represents a purely imaginary number (im * j).
    Imag { im: f32 },
    /// Represents a complex number with radius 1 and phase angle \(\phi \in (-\pi, \pi] \).
    Phase { ph: f32 },
    /// Represents a complex number in Cartesian coordinates (re + j * im).
    Ccs { re: f32, im: f32 },
    /// Represents a complex number in logarithmic form: \(\ln(z) = \text{re} + j \cdot \text{im} \),
    /// such that \(z = e^{\text{re} + j \cdot \text{im}} \).
    Ln { re: f32, im: f32 },
    /// Represents a complex number in polar coordinates: \( \text{rad} \cdot e^{j \cdot \text{ph}} \).
    PL { rad: f32, ph: f32 },
}

/// Enum representing error types for complex number operations.
#[derive(Debug, PartialEq)]
pub enum CpxError {
    /// Occurs when attempting to divide by zero.
    DivisionByZero,
    // Additional error types can be added as needed.
}

impl Hash for Cpx {
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

impl Neg for Cpx {
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
                ph: Cpx::branch_cut(ph + PI),
            },
            Cpx::Ccs { re, im } => Cpx::Ccs { re: -re, im: -im },
            Cpx::Ln { re, im } => Cpx::Ln {
                re,
                im: Cpx::branch_cut(im + PI),
            },
            Cpx::PL { rad, ph } => Cpx::PL {
                rad,
                ph: Cpx::branch_cut(ph + PI),
            },
        }
    }
}

impl PartialEq for Cpx {
    fn eq(&self, other: &Self) -> bool {
        match (self.regularize(), other.regularize()) {
            (Cpx::Zero {}, Cpx::Zero {}) => true,
            (Cpx::One {}, Cpx::One {}) => true,
            (Cpx::NegOne {}, Cpx::NegOne {}) => true,
            (Cpx::J {}, Cpx::J {}) => true,
            (Cpx::NegJ {}, Cpx::NegJ {}) => true,
            (Cpx::Real { re: re1 }, Cpx::Real { re: re2 }) => re1 == re2,
            (Cpx::Imag { im: im1 }, Cpx::Imag { im: im2 }) => im1 == im2,
            (Cpx::Phase { ph: ph1 }, Cpx::Phase { ph: ph2 }) => ph1 == ph2,
            (Cpx::Ccs { re: re1, im: im1 }, Cpx::Ccs { re: re2, im: im2 }) => {
                re1 == re2 && im1 == im2
            }
            (Cpx::Ln { re: re1, im: im1 }, Cpx::Ln { re: re2, im: im2 }) => {
                re1 == re2 && im1 == im2
            }
            (Cpx::PL { rad: rad1, ph: ph1 }, Cpx::PL { rad: rad2, ph: ph2 }) => {
                rad1 == rad2 && ph1 == ph2
            }
            _ => false,
        }
    }
}
impl Eq for Cpx {}

impl Cpx {
    /// Wraps a phase angle into the (-π, π] range.
    fn branch_cut(ph: f32) -> f32 {
        let wrapped = (ph + PI).rem_euclid(TAU) - PI;
        if wrapped <= -PI {
            wrapped + TAU
        } else {
            wrapped
        }
    }

    /// Returns a canonicalized form of the complex number, reducing floating-point errors and mapping
    /// close values to standard representations (e.g., 0, 1, -1, j, -j).
    pub fn regularize(self) -> Self {
        let threshold: f32 = 1e-6;
        match self {
            Cpx::Zero {} => Cpx::Zero {},
            Cpx::One {} => Cpx::One {},
            Cpx::NegOne {} => Cpx::NegOne {},
            Cpx::J {} => Cpx::J {},
            Cpx::NegJ {} => Cpx::NegJ {},
            Cpx::Real { re } => {
                if approx_eq(re, 0.0, threshold) {
                    Cpx::Zero {}
                } else if approx_eq(re, 1.0, threshold) {
                    Cpx::One {}
                } else if approx_eq(re, -1.0, threshold) {
                    Cpx::NegOne {}
                } else {
                    Cpx::Real { re }
                }
            }
            Cpx::Imag { im } => {
                if approx_eq(im, 0.0, threshold) {
                    Cpx::Zero {}
                } else if approx_eq(im, 1.0, threshold) {
                    Cpx::J {}
                } else if approx_eq(im, -1.0, threshold) {
                    Cpx::NegJ {}
                } else {
                    Cpx::Imag { im }
                }
            }
            Cpx::Phase { ph } => {
                let new_ph = Cpx::branch_cut(ph);
                if approx_eq(new_ph, 0.0, threshold) {
                    Cpx::One {}
                } else if approx_eq(new_ph, FRAC_PI_2, threshold) {
                    Cpx::J {}
                } else if approx_eq(new_ph, PI, threshold) {
                    Cpx::NegOne {}
                } else if approx_eq(new_ph, -FRAC_PI_2, threshold) {
                    Cpx::NegJ {}
                } else {
                    Cpx::Phase { ph: new_ph }
                }
            }
            Cpx::Ccs { re, im } => {
                let mag = re.hypot(im);
                let re_is_zero = approx_eq(re, 0.0, threshold);
                let im_is_zero = approx_eq(im, 0.0, threshold);

                if approx_eq(mag, 0.0, threshold) {
                    Cpx::Zero {}
                } else if approx_eq(re, 1.0, threshold) && im_is_zero {
                    Cpx::One {}
                } else if approx_eq(re, -1.0, threshold) && im_is_zero {
                    Cpx::NegOne {}
                } else if re_is_zero && approx_eq(im, 1.0, threshold) {
                    Cpx::J {}
                } else if re_is_zero && approx_eq(im, -1.0, threshold) {
                    Cpx::NegJ {}
                } else if im_is_zero {
                    Cpx::Real { re }
                } else if re_is_zero {
                    Cpx::Imag { im }
                } else if approx_eq(mag, 1.0, threshold) {
                    let phase = im.atan2(re);
                    Cpx::Phase { ph: phase }
                } else {
                    Cpx::Ccs { re, im }
                }
            }
            Cpx::Ln { re, im } => {
                let new_im = Cpx::branch_cut(im);
                let re_is_zero = approx_eq(re, 0.0, threshold);

                if re < threshold.ln() {
                    Cpx::Zero {}
                } else if re_is_zero {
                    match new_im {
                        x if approx_eq(x, 0.0, threshold) => Cpx::One {},
                        x if approx_eq(x, FRAC_PI_2, threshold) => Cpx::J {},
                        x if approx_eq(x, -FRAC_PI_2, threshold) => Cpx::NegJ {},
                        x if approx_eq(x, PI, threshold) => Cpx::NegOne {},
                        _ => Cpx::Phase { ph: new_im },
                    }
                } else if approx_eq(new_im, 0.0, threshold) {
                    Cpx::Real { re: re.exp() }
                } else if approx_eq(new_im, PI, threshold) {
                    Cpx::Real { re: -re.exp() }
                } else if approx_eq(new_im, FRAC_PI_2, threshold) {
                    Cpx::Imag { im: re.exp() }
                } else if approx_eq(new_im, -FRAC_PI_2, threshold) {
                    Cpx::Imag { im: -re.exp() }
                } else {
                    Cpx::Ln { re, im: new_im }
                }
            }
            Cpx::PL { rad, ph } => {
                let new_ph = Cpx::branch_cut(ph);
                let rad_is_zero = approx_eq(rad, 0.0, threshold);
                let rad_is_one = approx_eq(rad, 1.0, threshold);

                if rad_is_zero {
                    Cpx::Zero {}
                } else if rad_is_one && approx_eq(new_ph, 0.0, threshold) {
                    Cpx::One {}
                } else if rad_is_one && approx_eq(new_ph, FRAC_PI_2, threshold) {
                    Cpx::J {}
                } else if rad_is_one && approx_eq(new_ph, PI, threshold) {
                    Cpx::NegOne {}
                } else if rad_is_one && approx_eq(new_ph, -FRAC_PI_2, threshold) {
                    Cpx::NegJ {}
                } else if rad_is_one {
                    Cpx::Phase { ph: new_ph }
                } else if approx_eq(new_ph, 0.0, threshold) {
                    Cpx::Real { re: rad }
                } else if approx_eq(new_ph, PI, threshold) {
                    Cpx::Real { re: -rad }
                } else if approx_eq(new_ph, FRAC_PI_2, threshold) {
                    Cpx::Imag { im: rad }
                } else if approx_eq(new_ph, -FRAC_PI_2, threshold) {
                    Cpx::Imag { im: -rad }
                } else {
                    Cpx::PL { rad, ph: new_ph }
                }
            }
        }
    }

    /// Returns the complex conjugate of the complex number.
    pub fn conj(&self) -> Self {
        match self {
            Cpx::Zero {} => Cpx::Zero {},
            Cpx::One {} => Cpx::One {},
            Cpx::NegOne {} => Cpx::NegOne {},
            Cpx::J {} => Cpx::NegJ {},
            Cpx::NegJ {} => Cpx::J {},
            Cpx::Real { re } => Cpx::Real { re: *re },
            Cpx::Imag { im } => Cpx::Imag { im: -im },
            Cpx::Phase { ph } => Cpx::Phase {
                ph: Cpx::branch_cut(-ph),
            },
            Cpx::Ccs { re, im } => Cpx::Ccs { re: *re, im: -im },
            Cpx::Ln { re, im } => Cpx::Ln {
                re: *re,
                im: Cpx::branch_cut(-im),
            },
            Cpx::PL { rad, ph } => Cpx::PL {
                rad: *rad,
                ph: Cpx::branch_cut(-ph),
            },
        }
    }

    /// Returns the real part of the complex number as f32.
    pub fn re(&self) -> f32 {
        match self {
            Cpx::Zero {} => 0.0,
            Cpx::One {} => 1.0,
            Cpx::NegOne {} => -1.0,
            Cpx::J {} => 0.0,
            Cpx::NegJ {} => 0.0,
            Cpx::Real { re } => *re,
            Cpx::Imag { .. } => 0.0,
            Cpx::Phase { ph } => ph.cos(),
            Cpx::Ccs { re, .. } => *re,
            Cpx::Ln { re, im } => re.exp() * im.cos(),
            Cpx::PL { rad, ph } => rad * ph.cos(),
        }
    }

    /// Returns the imaginary part of the complex number as `f32`.
    pub fn im(&self) -> f32 {
        match self {
            Cpx::Zero {} => 0.0,
            Cpx::One {} => 0.0,
            Cpx::NegOne {} => 0.0,
            Cpx::J {} => 1.0,
            Cpx::NegJ {} => -1.0,
            Cpx::Real { .. } => 0.0,
            Cpx::Imag { im } => *im,
            Cpx::Phase { ph } => ph.sin(),
            Cpx::Ccs { im, .. } => *im,
            Cpx::Ln { re, im } => re.exp() * im.sin(),
            Cpx::PL { rad, ph } => rad * ph.sin(),
        }
    }

    /// Returns the radius (magnitude) of the complex number as `f32`.
    pub fn rad(&self) -> f32 {
        match self {
            Cpx::Zero {} => 0.0,
            Cpx::One {} => 1.0,
            Cpx::NegOne {} => 1.0,
            Cpx::J {} => 1.0,
            Cpx::NegJ {} => 1.0,
            Cpx::Real { re } => re.abs(),
            Cpx::Imag { im } => im.abs(),
            Cpx::Phase { .. } => 1.0,
            Cpx::Ccs { re, im } => re.hypot(*im),
            Cpx::Ln { re, .. } => re.exp(),
            Cpx::PL { rad, .. } => *rad,
        }
    }

    /// Returns the phase angle of the complex number in the interval (-π, π] as `f32`.
    pub fn ph(&self) -> f32 {
        match self {
            Cpx::Zero {} => 0.0,
            Cpx::One {} => 0.0,
            Cpx::NegOne {} => PI,
            Cpx::J {} => FRAC_PI_2,
            Cpx::NegJ {} => -FRAC_PI_2,
            Cpx::Real { .. } => 0.0,
            Cpx::Imag { .. } => FRAC_PI_2,
            Cpx::Phase { ph } => Cpx::branch_cut(*ph),
            Cpx::Ln { im, .. } => Cpx::branch_cut(*im),
            Cpx::PL { ph, .. } => Cpx::branch_cut(*ph),
            Cpx::Ccs { re, im } => im.atan2(*re),
        }
    }

    /// Returns a complex number representing the rotation (phase) of the input complex number.
    pub fn rot(&self) -> Cpx {
        Cpx::Phase { ph: self.ph() }
    }

    /// Returns the square root of the complex number.
    pub fn sqrt(&self) -> Self {
        match self {
            Cpx::Zero {} => Cpx::Zero {},
            Cpx::One {} => Cpx::One {},
            Cpx::NegOne {} => Cpx::J {},
            Cpx::J {} => Cpx::Phase { ph: FRAC_PI_4 },
            Cpx::NegJ {} => Cpx::Phase { ph: -FRAC_PI_4 },
            Cpx::Real { re } => {
                if *re >= 0.0 {
                    Cpx::Real { re: re.sqrt() }
                } else {
                    Cpx::Imag { im: (-re).sqrt() }
                }
            }
            Cpx::Imag { im } => {
                if *im >= 0.0 {
                    Cpx::PL {
                        rad: im.abs(),
                        ph: FRAC_PI_4,
                    }
                } else {
                    Cpx::PL {
                        rad: im.abs(),
                        ph: -FRAC_PI_4,
                    }
                }
            }
            Cpx::Phase { ph } => Cpx::Phase { ph: ph / 2.0 },
            Cpx::Ccs { re, im } => Cpx::Ccs { re: *re, im: -im },
            Cpx::Ln { re, im } => Cpx::Ln {
                re: re / 2.0,
                im: im / 2.0,
            },
            Cpx::PL { rad, ph } => Cpx::PL {
                rad: rad.sqrt(),
                ph: ph / 2.0,
            },
        }
    }

    /// Returns the exponential of the complex number.
    pub fn exp(&self) -> Self {
        match self {
            Cpx::Zero {} => Cpx::One {},
            Cpx::One {} => Cpx::Real { re: E },
            Cpx::NegOne {} => Cpx::Real { re: 1.0 / E },
            Cpx::J {} => Cpx::Phase { ph: 1.0 },
            Cpx::NegJ {} => Cpx::Phase { ph: -1.0 },
            Cpx::Real { re } => Cpx::Real { re: re.exp() },
            Cpx::Imag { im } => Cpx::Phase {
                ph: Cpx::branch_cut(*im),
            }
            .regularize(),
            Cpx::Phase { ph } => Cpx::Ln {
                re: ph.cos(),
                im: ph.sin(),
            }
            .regularize(),
            Cpx::Ccs { re, im } => Cpx::Ln {
                re: *re,
                im: Cpx::branch_cut(*im),
            },
            Cpx::Ln { re, im } => Cpx::Ln {
                re: re.exp() * im.cos(),
                im: Cpx::branch_cut(re.exp() * im.sin()),
            }
            .regularize(),
            Cpx::PL { rad, ph } => Cpx::Ln {
                re: rad * ph.cos(),
                im: Cpx::branch_cut(rad * ph.sin()),
            }
            .regularize(),
        }
    }

    /// Returns `true` if the complex number is zero, and `false` otherwise.
    pub fn is_zero(&self) -> bool {
        matches!(self.regularize(), Cpx::Zero {})
    }

    /// Returns the inverse (reciprocal) of the complex number.
    pub fn inv(&self) -> Result<Self, CpxError> {
        match self {
            Cpx::Zero {} => Err(CpxError::DivisionByZero),
            Cpx::One {} => Ok(Cpx::One {}),
            Cpx::NegOne {} => Ok(Cpx::NegOne {}),
            Cpx::J {} => Ok(Cpx::NegJ {}),
            Cpx::NegJ {} => Ok(Cpx::J {}),
            Cpx::Ccs { .. } => Ok(Cpx::PL {
                rad: 1.0 / self.rad(),
                ph: Cpx::branch_cut(-self.ph()),
            }),
            Cpx::Real { re } => Ok(Cpx::Real { re: 1.0 / re }),
            Cpx::Imag { im } => Ok(Cpx::Imag { im: -1.0 / im }),
            Cpx::Phase { ph } => Ok(Cpx::Phase {
                ph: Cpx::branch_cut(-ph),
            }),
            Cpx::Ln { re, im } => Ok(Cpx::Ln {
                re: -re,
                im: Cpx::branch_cut(-im),
            }),
            Cpx::PL { rad, ph } => Ok(Cpx::PL {
                rad: 1.0 / rad,
                ph: Cpx::branch_cut(-ph),
            }),
        }
    }
    /// Raises the complex number to an integer power using polar form.
    ///
    /// Computes `r^n * e^(i * n * θ)` by raising the magnitude to `n` and scaling the phase.
    /// Returns the regularized result.
    pub fn powi(self, n: i32) -> Self {
        if n == 0 {
            return ONE;
        }
        let new_rad = self.rad().powi(n);
        let new_ph = Cpx::branch_cut(self.ph() * (n as f32));
        Cpx::PL {
            rad: new_rad,
            ph: new_ph,
        }
    }
}

impl Add for Cpx {
    type Output = Cpx;
    fn add(self, other: Cpx) -> Cpx {
        match (self, other) {
            (Cpx::Zero {}, x) | (x, Cpx::Zero {}) => x, // Adding zero is a no-op
            (Cpx::One {}, Cpx::One {}) => Cpx::Real { re: 2.0 },
            (Cpx::One {}, Cpx::NegOne {}) | (Cpx::NegOne {}, Cpx::One {}) => Cpx::Zero {},
            (Cpx::One {}, Cpx::J {}) | (Cpx::J {}, Cpx::One {}) => Cpx::Ccs { re: 1.0, im: 1.0 },
            (Cpx::One {}, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::One {}) => {
                Cpx::Ccs { re: 1.0, im: -1.0 }
            }
            (Cpx::NegOne {}, Cpx::NegOne {}) => Cpx::Real { re: -2.0 },
            (Cpx::NegOne {}, Cpx::J {}) | (Cpx::J {}, Cpx::NegOne {}) => {
                Cpx::Ccs { re: -1.0, im: 1.0 }
            }
            (Cpx::NegOne {}, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::NegOne {}) => {
                Cpx::Ccs { re: -1.0, im: -1.0 }
            }
            (Cpx::J {}, Cpx::J {}) => Cpx::Imag { im: 2.0 },
            (Cpx::J {}, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::J {}) => Cpx::Zero {},
            (Cpx::NegJ {}, Cpx::NegJ {}) => Cpx::Imag { im: -2.0 },
            (Cpx::Real { re: r1 }, Cpx::Real { re: r2 }) => Cpx::Real { re: r1 + r2 },
            (Cpx::Real { re }, Cpx::One {}) | (Cpx::One {}, Cpx::Real { re }) => {
                Cpx::Real { re: re + 1.0 }
            }
            (Cpx::Real { re }, Cpx::NegOne {}) | (Cpx::NegOne {}, Cpx::Real { re }) => {
                Cpx::Real { re: re - 1.0 }
            }
            (Cpx::Real { re }, Cpx::J {}) | (Cpx::J {}, Cpx::Real { re }) => {
                Cpx::Ccs { re, im: 1.0 }
            }
            (Cpx::Real { re }, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::Real { re }) => {
                Cpx::Ccs { re, im: -1.0 }
            }
            (Cpx::Imag { im: i1 }, Cpx::Imag { im: i2 }) => Cpx::Imag { im: i1 + i2 },
            (Cpx::Real { re }, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::Real { re }) => {
                Cpx::Ccs { re, im }
            }
            (Cpx::One {}, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::One {}) => {
                Cpx::Ccs { re: 1.0, im }
            }
            (Cpx::NegOne {}, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::NegOne {}) => {
                Cpx::Ccs { re: -1.0, im }
            }
            (Cpx::J {}, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::J {}) => {
                Cpx::Imag { im: im + 1.0 }
            }
            (Cpx::NegJ {}, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::NegJ {}) => {
                Cpx::Imag { im: im - 1.0 }
            }
            (Cpx::Ccs { re: r1, im: i1 }, Cpx::Ccs { re: r2, im: i2 }) => Cpx::Ccs {
                re: r1 + r2,
                im: i1 + i2,
            },
            (Cpx::Ccs { re, im }, Cpx::One {}) | (Cpx::One {}, Cpx::Ccs { re, im }) => {
                Cpx::Ccs { re: re + 1.0, im }
            }
            (Cpx::Ccs { re, im }, Cpx::NegOne {}) | (Cpx::NegOne {}, Cpx::Ccs { re, im }) => {
                Cpx::Ccs { re: re - 1.0, im }
            }
            (Cpx::Ccs { re, im }, Cpx::J {}) | (Cpx::J {}, Cpx::Ccs { re, im }) => {
                Cpx::Ccs { re, im: im + 1.0 }
            }
            (Cpx::Ccs { re, im }, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::Ccs { re, im }) => {
                Cpx::Ccs { re, im: im - 1.0 }
            }
            (Cpx::Ccs { re, im }, Cpx::Real { re: r2 })
            | (Cpx::Real { re: r2 }, Cpx::Ccs { re, im }) => Cpx::Ccs { re: re + r2, im },
            (Cpx::Ccs { re, im }, Cpx::Imag { im: i2 })
            | (Cpx::Imag { im: i2 }, Cpx::Ccs { re, im }) => Cpx::Ccs { re, im: im + i2 },
            _ => Cpx::Ccs {
                re: self.re() + other.re(),
                im: self.im() + other.im(),
            },
        }
    }
}
impl AddAssign for Cpx {
    fn add_assign(&mut self, rhs: Self) {
        *self = *self + rhs;
    }
}
impl Sub for Cpx {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}
impl SubAssign for Cpx {
    fn sub_assign(&mut self, rhs: Self) {
        *self = *self - rhs;
    }
}
impl Mul for Cpx {
    type Output = Cpx;
    fn mul(self, other: Cpx) -> Cpx {
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
                Cpx::Phase { ph: ph + PI }.regularize()
            }
            (Cpx::Phase { ph }, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::Phase { ph }) => {
                Cpx::Phase { ph: ph - PI }.regularize()
            }
            (Cpx::Phase { ph }, Cpx::Real { re }) | (Cpx::Real { re }, Cpx::Phase { ph }) => {
                Cpx::PL { rad: re, ph }.regularize()
            }
            (Cpx::Phase { ph }, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::Phase { ph }) => {
                if im >= 0.0 {
                    Cpx::PL {
                        rad: im.abs(),
                        ph: ph + FRAC_PI_2,
                    }
                    .regularize()
                } else {
                    Cpx::PL {
                        rad: im.abs(),
                        ph: ph - FRAC_PI_2,
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
                im: im + FRAC_PI_2,
            }
            .regularize(),
            (Cpx::Ln { re, im }, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::Ln { re, im }) => Cpx::Ln {
                re,
                im: im - FRAC_PI_2,
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
                if im >= 0.0 {
                    Cpx::Ln {
                        re: re + i2.abs().ln(),
                        im: im + FRAC_PI_2,
                    }
                    .regularize()
                } else {
                    Cpx::Ln {
                        re: re + i2.abs().ln(),
                        im: im - FRAC_PI_2,
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
                ph: ph + FRAC_PI_2,
            }
            .regularize(),
            (Cpx::PL { rad, ph }, Cpx::NegJ {}) | (Cpx::NegJ {}, Cpx::PL { rad, ph }) => Cpx::PL {
                rad,
                ph: ph - FRAC_PI_2,
            }
            .regularize(),
            (Cpx::PL { rad, ph }, Cpx::Real { re }) | (Cpx::Real { re }, Cpx::PL { rad, ph }) => {
                Cpx::PL { rad: rad * re, ph }.regularize()
            }
            (Cpx::PL { rad, ph }, Cpx::Imag { im }) | (Cpx::Imag { im }, Cpx::PL { rad, ph }) => {
                if im >= 0.0 {
                    Cpx::PL {
                        rad: rad * im.abs(),
                        ph: ph + FRAC_PI_2,
                    }
                    .regularize()
                } else {
                    Cpx::PL {
                        rad: rad * im.abs(),
                        ph: ph - FRAC_PI_2,
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
impl MulAssign for Cpx {
    fn mul_assign(&mut self, rhs: Self) {
        *self = *self * rhs;
    }
}
impl Div for Cpx {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let new_rad = self.rad() / other.rad();
        let new_phase = Cpx::branch_cut(self.ph() - other.ph());
        Cpx::PL {
            rad: new_rad,
            ph: new_phase,
        }
        .regularize()
    }
}
impl DivAssign for Cpx {
    fn div_assign(&mut self, rhs: Self) {
        *self = *self / rhs;
    }
}
impl Add<f32> for Cpx {
    type Output = Self;
    fn add(self, other: f32) -> Self {
        self + Cpx::Real { re: other }
    }
}
impl AddAssign<f32> for Cpx {
    fn add_assign(&mut self, rhs: f32) {
        *self = *self + rhs;
    }
}
impl Add<Cpx> for f32 {
    type Output = Cpx;
    fn add(self, other: Cpx) -> Cpx {
        other + self
    }
}
impl Sub<f32> for Cpx {
    type Output = Self;
    fn sub(self, other: f32) -> Self {
        self + (-other)
    }
}
impl SubAssign<f32> for Cpx {
    fn sub_assign(&mut self, rhs: f32) {
        *self = *self - rhs;
    }
}
impl Sub<Cpx> for f32 {
    type Output = Cpx;
    fn sub(self, other: Cpx) -> Cpx {
        (-other) + self
    }
}
impl Mul<f32> for Cpx {
    type Output = Self;
    fn mul(self, other: f32) -> Self {
        self * Cpx::Real { re: other }
    }
}
impl MulAssign<f32> for Cpx {
    fn mul_assign(&mut self, rhs: f32) {
        *self = *self * rhs;
    }
}
impl Mul<Cpx> for f32 {
    type Output = Cpx;
    fn mul(self, other: Cpx) -> Cpx {
        other * self
    }
}
impl Div<f32> for Cpx {
    type Output = Self;
    fn div(self, other: f32) -> Self {
        self * (1.0 / other)
    }
}
impl DivAssign<f32> for Cpx {
    fn div_assign(&mut self, rhs: f32) {
        *self = *self / rhs;
    }
}
impl Div<Cpx> for f32 {
    type Output = Cpx;
    fn div(self, other: Cpx) -> Cpx {
        Cpx::Real { re: self } / other
    }
}

/// Approximately checks whether two `f32` numbers are equal within a given epsilon tolerance.
///
/// # Arguments
/// * `a` - First floating-point number.
/// * `b` - Second floating-point number.
/// * `eps` - Maximum allowed difference to consider `a` and `b` approximately equal.
///
/// # Returns
/// * `true` if the absolute difference between `a` and `b` is less than or equal to `eps`, otherwise `false`.
fn approx_eq(a: f32, b: f32, eps: f32) -> bool {
    (a - b).abs() <= eps
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_constants() {
        assert_eq!(ZERO, Cpx::Zero {});
        assert_eq!(ONE, Cpx::One {});
        assert_eq!(NEG_ONE, Cpx::NegOne {});
        assert_eq!(J, Cpx::J {});
        assert_eq!(NEG_J, Cpx::NegJ {});
        assert_eq!(INV_SQRT_2, Cpx::Real { re: 1.0 / SQRT_2 });
        assert_eq!(NEG_INV_SQRT_2, Cpx::Real { re: -1.0 / SQRT_2 });
        assert_eq!(J_INV_SQRT_2, Cpx::Imag { im: 1.0 / SQRT_2 });
        assert_eq!(NEG_J_INV_SQRT_2, Cpx::Imag { im: -1.0 / SQRT_2 });
        assert_eq!(SQRT_J, Cpx::Phase { ph: FRAC_PI_4 });
    }

    #[test]
    fn test_hash() {
        use std::collections::hash_map::DefaultHasher;
        fn calculate_hash<T: Hash>(t: &T) -> u64 {
            let mut s = DefaultHasher::new();
            t.hash(&mut s);
            s.finish()
        }

        assert_eq!(calculate_hash(&ZERO), calculate_hash(&Cpx::Zero {}));
        assert_eq!(calculate_hash(&ONE), calculate_hash(&Cpx::One {}));
        assert_eq!(calculate_hash(&Cpx::Real { re: 1.0 }), calculate_hash(&ONE));
        assert_ne!(calculate_hash(&Cpx::Real { re: 1.1 }), calculate_hash(&ONE));
    }

    #[test]
    fn test_neg() {
        assert_eq!(-ZERO, ZERO);
        assert_eq!(-ONE, NEG_ONE);
        assert_eq!(-NEG_ONE, ONE);
        assert_eq!(-J, NEG_J);
        assert_eq!(-NEG_J, J);
        assert_eq!(-Cpx::Real { re: 1.0 }, Cpx::Real { re: -1.0 });
        assert_eq!(-Cpx::Imag { im: 2.0 }, Cpx::Imag { im: -2.0 });
        assert_eq!(
            -Cpx::Phase { ph: FRAC_PI_4 },
            Cpx::Phase {
                ph: -3.0 * FRAC_PI_4
            }
        );
        assert_eq!(
            -Cpx::Ccs { re: 1.0, im: 2.0 },
            Cpx::Ccs { re: -1.0, im: -2.0 }
        );
        assert_eq!(
            -Cpx::Ln { re: 1.0, im: 2.0 },
            Cpx::Ln {
                re: 1.0,
                im: 2.0 + PI
            }
        );
        assert_eq!(
            -Cpx::PL { rad: 1.0, ph: 2.0 },
            Cpx::PL {
                rad: 1.0,
                ph: 2.0 + PI
            }
        );
    }

    #[test]
    fn test_eq() {
        assert_eq!(ZERO, ZERO);
        assert_eq!(ONE, ONE);
        assert_eq!(NEG_ONE, NEG_ONE);
        assert_eq!(J, J);
        assert_eq!(NEG_J, NEG_J);
        assert_eq!(Cpx::Real { re: 1.0 }, Cpx::Real { re: 1.0 });
        assert_eq!(Cpx::Imag { im: 2.0 }, Cpx::Imag { im: 2.0 });
        assert_eq!(Cpx::Phase { ph: FRAC_PI_4 }, Cpx::Phase { ph: FRAC_PI_4 });
        assert_eq!(Cpx::Ccs { re: 1.0, im: 2.0 }, Cpx::Ccs { re: 1.0, im: 2.0 });
        assert_eq!(Cpx::Ln { re: 1.0, im: 2.0 }, Cpx::Ln { re: 1.0, im: 2.0 });
        assert_eq!(Cpx::PL { rad: 1.0, ph: 2.0 }, Cpx::PL { rad: 1.0, ph: 2.0 });
        assert_eq!(Cpx::Real { re: 1.0 }, ONE);
        assert_eq!(Cpx::Imag { im: 1.0 }, J);
        assert_ne!(Cpx::Real { re: 1.0 }, Cpx::Imag { im: 1.0 });
        assert_ne!(Cpx::Real { re: 1.0 }, Cpx::Real { re: 1.1 });
    }

    #[test]
    fn test_regularize() {
        assert_eq!(Cpx::Real { re: 0.0 }, ZERO);
        assert_eq!(Cpx::Real { re: 1.0 }, ONE);
        assert_eq!(Cpx::Real { re: -1.0 }, NEG_ONE);
        assert_eq!(Cpx::Imag { im: 0.0 }, ZERO);
        assert_eq!(Cpx::Imag { im: 1.0 }, J);
        assert_eq!(Cpx::Imag { im: -1.0 }, NEG_J);
        assert_eq!(Cpx::Phase { ph: 0.0 }, ONE);
        assert_eq!(Cpx::Phase { ph: FRAC_PI_2 }, J);
        assert_eq!(Cpx::Phase { ph: PI }, NEG_ONE);
        assert_eq!(Cpx::Phase { ph: -FRAC_PI_2 }, NEG_J);
        assert_eq!(Cpx::Ccs { re: 0.0, im: 0.0 }, ZERO);
        assert_eq!(Cpx::Ccs { re: 1.0, im: 0.0 }, ONE);
        assert_eq!(Cpx::Ccs { re: 0.0, im: 1.0 }, J);
        assert_eq!(Cpx::Ln { re: 0.0, im: 0.0 }, ONE);
        assert_eq!(Cpx::PL { rad: 0.0, ph: 0.0 }, ZERO);
        assert_eq!(Cpx::PL { rad: 1.0, ph: 0.0 }, ONE);
        assert_eq!(
            Cpx::PL {
                rad: 1.0,
                ph: FRAC_PI_2
            },
            J
        );
    }
    #[test]
    fn test_conj() {
        assert_eq!(ZERO.conj(), ZERO);
        assert_eq!(ONE.conj(), ONE);
        assert_eq!(J.conj(), NEG_J);
        assert_eq!(Cpx::Real { re: 1.0 }.conj(), Cpx::Real { re: 1.0 });
        assert_eq!(Cpx::Imag { im: 2.0 }.conj(), Cpx::Imag { im: -2.0 });
        assert_eq!(
            Cpx::Phase { ph: FRAC_PI_4 }.conj(),
            Cpx::Phase { ph: -FRAC_PI_4 }
        );
        assert_eq!(
            Cpx::Ccs { re: 1.0, im: 2.0 }.conj(),
            Cpx::Ccs { re: 1.0, im: -2.0 }
        );
        assert_eq!(
            Cpx::Ln { re: 1.0, im: 2.0 }.conj(),
            Cpx::Ln { re: 1.0, im: -2.0 }
        );
        assert_eq!(
            Cpx::PL { rad: 1.0, ph: 2.0 }.conj(),
            Cpx::PL { rad: 1.0, ph: -2.0 }
        );
    }

    #[test]
    fn test_re_im_rad_ph() {
        assert_eq!(ONE.re(), 1.0);
        assert_eq!(J.im(), 1.0);
        assert_eq!(Cpx::Ccs { re: 3.0, im: 4.0 }.rad(), 5.0);
        assert_eq!(J.ph(), FRAC_PI_2);
    }
    #[test]
    fn test_rot() {
        assert_eq!(
            Cpx::Ccs { re: 3.0, im: 4.0 }.rot(),
            Cpx::Phase {
                ph: 4.0f32.atan2(3.0)
            }
        );
    }

    #[test]
    fn test_sqrt() {
        assert_eq!(ONE.sqrt(), ONE);
        assert_eq!(NEG_ONE.sqrt(), J);
        assert_eq!(J.sqrt(), SQRT_J);
        assert_eq!(Cpx::Real { re: 4.0 }.sqrt(), Cpx::Real { re: 2.0 });
        assert_eq!(
            Cpx::Imag { im: 4.0 }.sqrt(),
            Cpx::PL {
                rad: 4.0,
                ph: FRAC_PI_4
            }
        );
        assert_eq!(Cpx::Phase { ph: PI }.sqrt(), Cpx::Phase { ph: FRAC_PI_2 });
        assert_eq!(
            Cpx::PL { rad: 4.0, ph: PI }.sqrt(),
            Cpx::PL {
                rad: 2.0,
                ph: FRAC_PI_2
            }
        );
    }
    #[test]
    fn test_exp() {
        assert_eq!(ZERO.exp(), ONE);
        assert_eq!(ONE.exp(), Cpx::Real { re: E });
        assert_eq!(J.exp(), Cpx::Phase { ph: 1.0 });
        //assert_eq!(Cpx::Real { re: 2.0 }.exp(), Cpx::Real { re: E * E });
        // Use a tolerance for floating-point comparisons
        let expected = Cpx::Real { re: E * E };
        let actual = Cpx::Real { re: 2.0 }.exp();
        let tolerance = 1e-6; // Adjust tolerance as needed

        match (actual, expected) {
            (Cpx::Real { re: a }, Cpx::Real { re: b }) => {
                assert!((a - b).abs() < tolerance);
            }
            _ => assert!(false, "Types don't match or not Real"),
        }
        assert_eq!(Cpx::Imag { im: PI }.exp(), NEG_ONE);
        assert_eq!(Cpx::Phase { ph: PI }.exp(), Cpx::Ln { re: -1.0, im: 0.0 });
    }

    #[test]
    fn test_is_zero() {
        assert!(ZERO.is_zero());
        assert!(!ONE.is_zero());
    }

    #[test]
    fn test_inv() {
        assert_eq!(ZERO.inv(), Err(CpxError::DivisionByZero));
        assert_eq!(ONE.inv(), Ok(ONE));
        assert_eq!(J.inv(), Ok(NEG_J));
        assert_eq!(Cpx::Real { re: 2.0 }.inv(), Ok(Cpx::Real { re: 0.5 }));
        assert_eq!(Cpx::Imag { im: 2.0 }.inv(), Ok(Cpx::Imag { im: -0.5 }));
        assert_eq!(
            Cpx::Ccs { re: 3.0, im: 4.0 }.inv(),
            Ok(Cpx::PL {
                rad: 0.2,
                ph: -4.0f32.atan2(3.0)
            })
        );
    }

    #[test]
    fn test_add() {
        assert_eq!(ONE + ONE, Cpx::Real { re: 2.0 });
        assert_eq!(ONE + J, Cpx::Ccs { re: 1.0, im: 1.0 });
        assert_eq!(
            Cpx::Real { re: 2.0 } + Cpx::Imag { im: 3.0 },
            Cpx::Ccs { re: 2.0, im: 3.0 }
        );
        assert_eq!(
            Cpx::Ccs { re: 1.0, im: 2.0 } + Cpx::Ccs { re: 3.0, im: 4.0 },
            Cpx::Ccs { re: 4.0, im: 6.0 }
        );
    }

    #[test]
    fn test_add_assign() {
        let mut c = ONE;
        c += ONE;
        assert_eq!(c, Cpx::Real { re: 2.0 });
    }
    #[test]
    fn test_sub() {
        assert_eq!(ONE - ONE, ZERO);
        assert_eq!(ONE - J, Cpx::Ccs { re: 1.0, im: -1.0 });
        assert_eq!(
            Cpx::Real { re: 2.0 } - Cpx::Imag { im: 3.0 },
            Cpx::Ccs { re: 2.0, im: -3.0 }
        );
        assert_eq!(
            Cpx::Ccs { re: 1.0, im: 2.0 } - Cpx::Ccs { re: 3.0, im: 4.0 },
            Cpx::Ccs { re: -2.0, im: -2.0 }
        );
    }

    #[test]
    fn test_sub_assign() {
        let mut c = ONE;
        c -= ONE;
        assert_eq!(c, ZERO);
    }

    #[test]
    fn test_mul() {
        assert_eq!(ONE * ONE, ONE);
        assert_eq!(J * J, NEG_ONE);
        assert_eq!(
            Cpx::Real { re: 2.0 } * Cpx::Imag { im: 3.0 },
            Cpx::Imag { im: 6.0 }
        );
        assert_eq!(
            Cpx::Phase { ph: FRAC_PI_4 } * Cpx::Phase { ph: FRAC_PI_4 },
            J
        );
        assert_eq!(
            Cpx::PL {
                rad: 2.0,
                ph: FRAC_PI_4
            } * Cpx::PL {
                rad: 3.0,
                ph: FRAC_PI_4
            },
            Cpx::PL {
                rad: 6.0,
                ph: FRAC_PI_2
            }
        );
    }

    #[test]
    fn test_mul_assign() {
        let mut c = ONE;
        c *= ONE;
        assert_eq!(c, ONE);
    }

    #[test]
    fn test_div() {
        assert_eq!(ONE / ONE, ONE);
        assert_eq!(J / J, ONE);
        assert_eq!(
            Cpx::Real { re: 4.0 } / Cpx::Real { re: 2.0 },
            Cpx::Real { re: 2.0 }
        );
        assert_eq!(
            Cpx::Imag { im: 6.0 } / Cpx::Real { re: 2.0 },
            Cpx::Imag { im: 3.0 }
        );
        assert_eq!(
            Cpx::PL {
                rad: 6.0,
                ph: FRAC_PI_2
            } / Cpx::PL {
                rad: 2.0,
                ph: FRAC_PI_4
            },
            Cpx::PL {
                rad: 3.0,
                ph: FRAC_PI_4
            }
        );
    }

    #[test]
    fn test_div_assign() {
        let mut c = ONE;
        c /= ONE;
        assert_eq!(c, ONE);
    }

    #[test]
    fn test_add_f32() {
        assert_eq!(ONE + 1.0, Cpx::Real { re: 2.0 });
        assert_eq!(1.0 + ONE, Cpx::Real { re: 2.0 });
    }
    #[test]
    fn test_add_f32_assign() {
        let mut c = ONE;
        c += 1.0;
        assert_eq!(c, Cpx::Real { re: 2.0 });
    }
    #[test]
    fn test_sub_f32() {
        assert_eq!(ONE - 1.0, ZERO);
        assert_eq!(2.0 - ONE, Cpx::Real { re: 1.0 });
    }
    #[test]
    fn test_sub_f32_assign() {
        let mut c = ONE;
        c -= 1.0;
        assert_eq!(c, ZERO);
    }
    #[test]
    fn test_mul_f32() {
        assert_eq!(ONE * 2.0, Cpx::Real { re: 2.0 });
        assert_eq!(2.0 * ONE, Cpx::Real { re: 2.0 });
    }
    #[test]
    fn test_mul_f32_assign() {
        let mut c = ONE;
        c *= 2.0;
        assert_eq!(c, Cpx::Real { re: 2.0 });
    }

    #[test]
    fn test_div_f32() {
        assert_eq!(Cpx::Real { re: 4.0 } / 2.0, Cpx::Real { re: 2.0 });
        assert_eq!(4.0 / Cpx::Real { re: 2.0 }, Cpx::Real { re: 2.0 });
    }
    #[test]
    fn test_div_f32_assign() {
        let mut c = Cpx::Real { re: 4.0 };
        c /= 2.0;
        assert_eq!(c, Cpx::Real { re: 2.0 });
    }
}
