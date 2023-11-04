# Enumerative Sphere Shaping Library

This project implements the [enumerative sphere shaping (ESS)](https://www.sps.tue.nl/ictlab/pdf/publications/frans/WillemsWuijts1993VTCBenelux.pdf) (implementation follows [this paper](https://doi.org/10.1109/TWC.2019.2951139)) and [optimum ESS (OESS)](https://doi.org/10.1109/JLT.2022.3201901) algorithms in Rust.
For ease of use it also contains the subproject `pyrsess` which provides Python bindings for the Rust code.

- [Documentation](https://docs.rs/rsess/)
- [Source](https://github.com/kit-cel/rsess/tree/main)
- [PyRSESS](https://github.com/kit-cel/rsess/tree/main/pyrsess)
- [Examples](https://github.com/kit-cel/rsess_examples)

## Installation

The Rust code can be compiled and run with `cargo run`.
To use the Python bindings refer to the README in the `pyrsess` subfolder.
An optimized build can be created using `cargo build --release`.

As RSESS uses the rug crate which in turn uses GMP, which can not trivially be installed on Windows, compiling on Windows might not work.

The documentation can be compiled with `cargo doc`.
It can then be found as `.html` files in `./target/doc/ess/`.
The entry point is `./target/doc/rsess/index.html`.
By default the documentation focuses on the public interface of `rsess`.
If the reader is interested in the inner workings of `rsess`, running `cargo doc --document-private-items` may yield additional insights.

## Testing

Some test are located in `src/tests.rs`, these can be run with `cargo test`.

## Rust Code Overview

* File `src/lib.rs`
	* Defines the trait `DistributionMatcher`
		* Fixes common methods for `ESS` and `OESS`
	* Defines the structs `ESS` and `OESS` which implement `DistributionMatcher`
		* User interfaces for the ESS / OESS algorithms
			* Offer the `encode` / `decode` functions
		* Offer utilities like calculating the amplitude distribution or average energy
	* Defines the struct `ASK`
		* Represents the amplitude shift keying modulation scheme
* File `src/trellis.rs`
	* Defines the struct `Trellis`
		* Represents the trellis used internally by the ESS and OESS algorithms
		* Implements (forward and reverse) trellis construction
		* Implements indexing algorithms used in `encode` / `decode`
* File `src/iterators.rs`
	* Defines an iterator `Amplitudes`
		* Iterates trough amplitude values in a given ASK with some extra constraints
	* Defines an iterator `Energies`
		* Iterates through energy levels in a trellis
* File `src/tests.rs`
	* Defines tests for the remaining code
	* Roughly divided into three sections: Trellis, ESS and OESS

## Citing this work

If you use this library in your research you can reference the following publication.

Frederik Ritter, Andrej Rode, and Laurent Schmalen. "Introducing RSESS: An open source enumerative sphere shaping implementation coded in Rust." Proceedings of the GNU Radio Conference. Vol. 8. No. 1. 2023.
```
@inproceedings{ritter2023introducing,
  title={Introducing RSESS: An open source enumerative sphere shaping implementation coded in Rust},
  author={Ritter, Frederik and Rode, Andrej and Schmalen, Laurent},
  booktitle={Proceedings of the GNU Radio Conference},
  volume={8},
  number={1},
  year={2023}
}
```
