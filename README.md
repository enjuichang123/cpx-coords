# cpx-coords

To power the rust-quantum project, this library provides a robust `Cpx` type, specifically optimized for the intensive complex number multiplications required for quantum gate applications and tensor product operations. It supports various coordinate representations to maximize performance.

## Usage

Add the following to your `Cargo.toml` (to the latest version):

```toml
[dependencies]
cpx-coords = "0.1.5"
```

# Features

- Multiple Coordinate Representations: To optimize performance, particularly for multiplication, this library supports logarithmic coordinates. Multiplication in logarithmic form reduces to addition, which is faster than Cartesian multiplication. Additionally, common-case optimizations further improve time complexity.
- Precision: Uses f32 or f64 for floating-point precision.
- It can be a key of BTreeMap, HashMap and so on. Sometimes, we want to represent a superposition state as a map of (separable_state,cpx_value). Here, a separable_state use CpxKey can be treated as a key of the desired map.

# Future direction of development

A complex-number analog of IEEE standard of floating numbers is designing. Consider there is one sign bit, a few exp bits, and a few significant bits for a real floating number. We want to interpret the sign bit as a special case of phase bits, while the remaining exp bits and significant bits describe the radius.
In this interpretation, we may save a few bits for a complex number in polar coordinate.

## License

This project is licensed under either of

- [MIT license](LICENSE-MIT)
- [Apache License, Version 2.0](LICENSE-APACHE)

at your option.


### Contributions

Unless you explicitly state otherwise, any contribution intentionally submitted for inclusion in this project by you shall be licensed as above, without any additional terms or conditions.
