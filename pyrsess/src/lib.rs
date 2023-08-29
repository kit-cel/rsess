use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;
use numpy::{IntoPyArray, PyArray, PyArray1, PyArray2};
use duplicate::duplicate_item;

use rsess::ASK;
use rsess::DistributionMatcher;

use rug::Integer;

/// Encoder / decoder using the enumerative sphere shaping algorithm
///
/// The implementation follows <https://doi.org/10.1109/TWC.2019.2951139>.
#[pyclass]
pub struct ESS {
    ess: rsess::ESS,
}

/// Encoder / decoder using the optimum enumerative sphere shaping algorithm
///
/// Unlike ESS, OESS needs its `e_max` parameter to be 'optimal'. This is the case, if there is no
/// lower `e_max` parameter which leads to a shaper capable of encoding the same number of bits.
/// An optimal `e_max` parameter can be found using `OESS.optimal_e_max()`.
///
/// The implementation follows <https://doi.org/10.1109/JLT.2022.3201901>.
#[pyclass]
pub struct OESS {
    oess: rsess::OESS,
}

#[pymethods]
impl OESS {
    /// optimal_e_max(e_max, n_max, ask, /)
    /// --
    ///
    /// Returns the next lower optimal `e_max` value
    ///
    /// An ESS encoder with the provided arguments can encode a certain number of bits. This
    /// function returns the lowest `e_max` parameter so that an ESS encoder with the new
    /// `e_max` stil encodes the same number of bits.
    #[staticmethod]
    pub fn optimal_e_max(e_max: usize, n_max: usize, ask_num: usize) -> PyResult<usize> {
        Ok(rsess::OESS::optimal_e_max(e_max, n_max, &ASK::new(ask_num)))
    }
}

#[duplicate_item(
    shaper   shaper_class;
    [ess]    [ESS];
    [oess]   [OESS];
)]
#[pymethods]
impl shaper_class {

    /// new(e_max, n_max, ask, /)
    /// --
    ///
    /// Returns an new instance with the given parameters.
    ///
    /// Args:
    ///     `e_max` - positive integer, maximum sequence energy
    ///     `n_max` - positive integer, sequence length
    ///     `ask`   - positive integer power of two, e.g. 8 for 8-ASK
    #[new]
    pub fn new(
        e_max: usize,
        n_max: usize,
        ask: usize,
    ) -> PyResult<Self>{
        let shaper = rsess::shaper_class::new(e_max, n_max, ASK::new(ask));
        Ok(shaper_class{shaper})
    }
    /// Returns the number of amplitudes in an amplitude sequence
    pub fn n_max(&self) -> PyResult<usize> {
        Ok(self.shaper.n_max())
    }
    /// Returns the maximum allowed energy for an amplitude sequence
    pub fn e_max(&self) -> PyResult<usize> {
        Ok(self.shaper.e_max())
    }
    /// Returns the number of the used ASK, e.g. 8 if 8-ASK is used
    pub fn ask_num(&self) -> PyResult<usize> {
        Ok(self.shaper.ask_num())
    }

    /// encode(index_bits, /)
    /// --
    ///
    /// Returns the amplitude sequence for the given bits as a numpy array
    ///
    /// The values in `index_bits` should be either `1` or `0`. (Currently other values
    /// are possible but this may change.)
    ///
    /// This function raises an exception if `index_bits` is invalid.
    ///
    /// Args:
    ///     `index_bits` - numpy array or list
    pub fn encode<'py>(&self, py: Python<'py>, index_bits: Vec<u128>) -> PyResult<&'py PyArray1<usize>> {

        // convert vec of index bits to Integer
        let index = index_bits.into_iter().fold(Integer::new(), |integer, bit| {
            (integer << 1) + bit
        });

        let sequence = self.shaper.encode(&index);
        match sequence {
            Ok(amplitudes) => Ok(amplitudes.into_pyarray(py)),
            Err(_) => Err(PyValueError::new_err("Invalid index_bits!")),
        }
    }

    /// multi_encode(multi_index_bits, /)
    /// --
    ///
    /// Returns the amplitude sequences for multiple given bit strings as a 2D numpy array
    ///
    /// The values in `multi_index_bits` should be either `1` or `0`.
    ///
    /// This function raises an exception if one of the index bit strings in `multi_index_bits` is invalid.
    ///
    /// Args:
    ///     `index_bits` - 2D numpy array
    pub fn multi_encode<'py>(&self, py: Python<'py>, multi_index_bits: Vec<Vec<u32>>) -> PyResult<&'py PyArray2<usize>> {
        let mut sequences: Vec<Vec<usize>> = Vec::with_capacity(multi_index_bits.len());
        for index_bits in multi_index_bits {
            // convert vec of index bits to Integer
            let index = index_bits.into_iter().fold(Integer::new(), |integer, bit| {
                (integer << 1) + bit
            });

            let sequence = self.shaper.encode(&index);
            if sequence.is_err() {
                return Err(PyValueError::new_err("Invalid index_bits!"))
            } else {
                sequences.push(sequence.expect("Can't be err in this if branch."))
            }
        }
        let arr = PyArray::from_vec2(py, &sequences).expect("Should be valid ndarray");
        Ok(arr)
    }

    /// decode(sequence, /)
    /// --
    ///
    /// Returns the index corresponding to the provided amplitude sequence as a numpy
    /// array of `1`s and `0`s
    ///
    /// Raises an exception if `sequence` is invalid.
    ///
    /// Args:
    ///     `sequence` - numpy array or list
    pub fn decode<'py>(&self, py: Python<'py>, sequence: Vec<usize>) -> PyResult<&'py PyArray1<u32>> {
        // decodes a sequence of amplitudes to a bit array with same length specified by `get_num_bits`

        let index = self.shaper.decode(&sequence);
        let index = match index {
            Ok(idx) => Ok(idx),
            Err(_) => Err(PyValueError::new_err("Invalid amplitude sequence!")),
        };
        let index = index?;

        // convert index to numpy array
        let len = self.shaper.num_data_bits() as usize;
        let mut bits_vec = vec![0; len];
        let mut mask = Integer::from(1);
        let zero = Integer::from(0);
        for i in (0..len).rev() {
            let masked_index: Integer = (&mask & &index).into();
            if &masked_index != &zero {
                bits_vec[i] = 1;
            }
            mask = mask << 1;
        }
        Ok(bits_vec.into_pyarray(py))
    }

    /// multi_decode(sequences, /)
    /// --
    ///
    /// Returns the indexes corresponding to the provided amplitude sequence as a 2D numpy
    /// array of `1`s and `0`s
    ///
    /// Raises an exception if any amplitude sequence in `sequences` is invalid.
    ///
    /// Args:
    ///     `sequences` - 2D numpy array
    pub fn multi_decode<'py>(&self, py: Python<'py>, sequences: Vec<Vec<usize>>) -> PyResult<&'py PyArray2<u32>> {
        let mut bit_vectors = Vec::with_capacity(sequences.len());

        for sequence in sequences {
            let index = self.shaper.decode(&sequence);
            let index = match index {
                Ok(idx) => idx,
                Err(_) => Integer::from(0), // Ignore the error
            };

            // convert index to numpy array
            let len = self.shaper.num_data_bits() as usize;
            let mut bits_vec = vec![0; len];
            let mut mask = Integer::from(1);
            let zero = Integer::from(0);
            for i in (0..len).rev() {
                let masked_index: Integer = (&mask & &index).into();
                if &masked_index != &zero {
                    bits_vec[i] = 1;
                }
                mask = mask << 1;
            }
            bit_vectors.push(bits_vec);
        }
        Ok(PyArray::from_vec2(py, &bit_vectors).unwrap())
    }

    /// Returns the number of bits encoded per amplitude sequence
    pub fn num_data_bits(&self) -> PyResult<u32> {
        Ok(self.shaper.num_data_bits())
    }
    /// Returns the probabilities of the amplitude values
    ///
    /// The probabilities are returned as an array with the lowest index corresponding to the
    /// lowest amplitude.
    pub fn amplitude_distribution<'py>(&self, py: Python<'py>) -> PyResult<&'py PyArray1<f32>> {
        Ok(self.shaper.amplitude_distribution().into_pyarray(py))
    }
    /// Returns the  probabilities of the amplitude sequence energies
    ///
    /// The probabilities are returned as an array. The corresponding energy values can be
    /// calculated with `n_max + 8i` (`n_max` is the number of amplitudes per sequence and `i`
    /// is the index in the returned array).
    pub fn energy_distribution<'py>(&self, py: Python<'py>) -> PyResult<&'py PyArray1<f32>> {
        Ok(self.shaper.energy_distribution().into_pyarray(py))
    }
    /// Returns the average energy of amplitude sequences
    pub fn average_energy(&self) -> PyResult<f32> {
        Ok(self.shaper.average_energy())
    }
    /// Returns the number of used amplitude sequences as a string
    ///
    /// This will always be a power of two to have a fixed number of bits per amplitude sequence.
    pub fn num_sequences(&self) -> PyResult<String>{
        Ok(self.shaper.num_sequences().to_string_radix(10))
    }
    /// Returns the maximum number of possible amplitude sequences as a string
    ///
    /// WARNING: Effect of limiting the used indexes to a power of two is not regarded!!!
    /// If it should be, use `num_sequences` instead.
    pub fn num_sequences_possible(&self) -> PyResult<String>{
        Ok(self.shaper.num_sequences_possible().to_string_radix(10))
    }
}



/// Python distribution matcher module implemented in Rust.
#[pymodule]
fn pyrsess(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<ESS>()?;
    m.add_class::<OESS>()?;
    Ok(())
}
