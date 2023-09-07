# Python Bindings for RSESS

Python bindings for the ESS and OESS algorithms implemented in Rust.
The bindings are created using [PyO3](https://github.com/PyO3/pyo3).

- [Documentation](https://kit-cel.github.io/rsess/pyrsess.html)
- [Source](https://github.com/kit-cel/rsess)
- [Examples](https://github.com/kit-cel/rsess_examples)

## Installation

There may currently be some issues installing PyRSESS on Windows as RSESS uses GMP which can not easily be built on Windows.

### Using PIP

Type `pip install pyrsess` into your favorite command line.

### From Source

1. Make sure that Rust and its package manager `cargo` are installed
2. Clone this repository with `git clone https://github.com/kit-cel/rsess.git`
3. Create a virtual python environment in a folder of your choice (e.g. `python -m venv $VENV_NAME`)
4. Activate the virtual environment (e.g. `cd $VENV_NAME; source bin/activate` if you are using Bash)
5. Install the `pyrsess` package with `pip`: `pip install $YOUR_PATH_TO/rsess/pyrsess`
	- If this fails, your `pip` may be to old. Try `pip install --upgrade pip`

## Development

Building can be done according to: https://pyo3.rs/v0.17.3/getting_started.html

TLDR: `pip install maturin; maturin develop`

## Project content

- `src/lib.rs`: PyO3 Rust to python bindings
- `pyess.pyi`: Python function type hints
