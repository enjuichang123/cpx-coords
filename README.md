# cpx-coords

To power the rust-quantum project, this library provides a robust `Cpx` type, specifically optimized for the intensive complex number multiplications required for quantum gate applications and tensor product operations. It supports various coordinate representations to maximize performance.

## Usage

Add this to your `Cargo.toml`:

```toml
[dependencies]
cpx-coords = "0.1.0" # Replace with the latest version
```

# Features

- Multiple Coordinate Representations: To optimize performance, particularly for multiplication, this library supports logarithmic coordinates. Multiplication in logarithmic form reduces to addition, which is faster than Cartesian multiplication. Additionally, common-case optimizations further improve time complexity.
- Regularization: Implements a `regularize()` method to normalize complex numbers, ensuring consistent representations and handling numerical edge cases.
- Precision: Uses f32 for floating-point precision.
- Comprehensive Operations: Implements standard arithmetic operations (addition, subtraction, multiplication, division), negation, conjugation, square root, exponentiation, and more.
- Constants: Provides a suite of predefined constants for common complex numbers.

# Coordinate Representations

This library uses the following coordinate representations for complex numbers:

-   `Zero {}`: Represents the additive identity (0).
-   `One {}`: Represents the multiplicative identity (1).
-   `NegOne {}`: Represents the negation of `One` (-1).
-   `J {}`: Represents the imaginary unit (j), a square root of -1.
-   `NegJ {}`: Represents the negation of `J` (-j).
-   `Real { re: f32 }`: Represents a real number (re).
-   `Imag { im: f32 }`: Represents a purely imaginary number (im * j).
-   `Phase { ph: f32 }`: Represents a complex number with a radius of 1 and a phase angle (ph) in the interval (-π, π].
-   `Ccs { re: f32, im: f32 }`: Represents a complex number in Cartesian coordinates (re + j * im).
-   `Ln { re: f32, im: f32 }`: Represents a complex number in logarithmic form (exp(re + j * im)).
-   `PL { rad: f32, ph: f32 }`: Represents a complex number in polar coordinates (rad * exp(j * ph)).

## License

This project is licensed under either of

- [MIT license](LICENSE-MIT.txt)
- [Apache License, Version 2.0](LICENSE-APACHE.txt)

at your option.


### Contributions

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this project by you shall be licensed as above, without any additional terms or conditions.
