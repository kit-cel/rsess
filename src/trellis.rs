use rug::Integer;
use rug::Rational;
use std::cmp;

use crate::iterators::{Amplitudes, Energies};
use crate::ASK;

/// [Trellis] is a bounded energy trellis with some utility functions
///
/// For details on what a bounded energy trellis is see section III-B in
/// <https://doi.org/10.1109/TWC.2019.2951139>.
pub struct Trellis {
    pub ask: ASK,
    pub e_max: usize,
    pub n_max: usize,
    data: Vec<Vec<Integer>>,
}

impl Trellis {
    /// Create a new [Trellis] instance
    ///
    /// The energy trellis is calculated internally and the returned [Trellis] can be used without
    /// further configuration.
    pub fn new(e_max: usize, n_max: usize, ask: ASK) -> Trellis {
        let data = Trellis::calculate_trellis(e_max, n_max, &ask, false);
        Trellis {
            ask,
            e_max,
            n_max,
            data,
        }
    }

    /// Create a new [Trellis] instance representing a partial trellis
    ///
    /// A partial trellis is an energy trellis, which only contains sequences which have the
    /// maximum allowed energy. For details see section IV-B in <https://doi.org/10.1109/JLT.2022.3201901>.
    ///
    /// The partial energy trellis is calculated internally and the returned [Trellis] can be used without
    /// further configuration.
    pub fn new_partial(e_max: usize, n_max: usize, ask: ASK) -> Trellis {
        let data = Trellis::calculate_trellis(e_max, n_max, &ask, true);
        Trellis {
            ask,
            e_max,
            n_max,
            data,
        }
    }

    /// Constructs the (optionally partial) energy trellis from passed parameters
    ///
    /// This is a low level function that [Trellis::new()] and [Trellis::new_partial()] use internally,
    /// it is not advised to use it directly.
    ///
    /// The trellis calculation is based on section III-B in <https://doi.org/10.1109/TWC.2019.2951139>.
    ///
    /// Optionally a partial trellis can be calculated, this is based on section IV-B
    /// in <https://doi.org/10.1109/JLT.2022.3201901>.
    fn calculate_trellis(
        e_max: usize,
        n_max: usize,
        ask: &ASK,
        partial: bool,
    ) -> Vec<Vec<Integer>> {
        let max_energy_step = (e_max - n_max) / 8;
        let mut data = vec![vec![Integer::from(0); 1 + max_energy_step]; 1 + n_max];

        for n in (0..n_max + 1).rev() {
            for e in Energies::for_trellis_level(n, n_max, e_max) {
                if n == n_max {
                    if partial {
                        // number of possible sequences for all end nodes is 0,
                        // only for the maximum energy end node it is 1
                        if e + 8 > e_max {
                            Trellis::set(&mut data, n, e, Integer::from(1));
                        } else {
                            Trellis::set(&mut data, n, e, Integer::from(0));
                        }
                    } else {
                        // number of possible sequences for all end nodes is 1
                        Trellis::set(&mut data, n, e, Integer::from(1));
                    }
                } else {
                    // number of possible paths for a node is the sum of the number
                    // of possible sequences of all successor nodes
                    let mut num_possible_sequences = Integer::from(0);
                    for u in Amplitudes::possible_from_trellis_node(n, e, n_max, e_max, ask) {
                        let successor_energy = u * u + e;
                        num_possible_sequences += Trellis::get(&data, n + 1, successor_energy);
                    }
                    Trellis::set(&mut data, n, e, num_possible_sequences);
                }
                // debugging output
                //println!("n: {}, e: {}, value: {}", n, node_energy, Trellis::get(&data, n, node_energy));
            }
        }

        data
    }

    /// Calculate the reverse trellis of a fully utilized trellis
    ///
    /// A number at location (`n`, `e`) in the reverse trellis is the number of length `n` sequences with
    /// the energy `e`. Most importantly a number at (`n_max`, `e`) signifies the number of sequences with
    /// energy `e` in an energy trellis with the same parameters as this function.
    ///
    /// The returned reverse trellis is in the shape of a 2D [Vec]. To access values via their `n`
    /// and `e` indexes use the [Trellis::get()] function.
    ///
    /// This is a simpler case of the algorithm in [Trellis::reverse_trellis()] because all
    /// sequences in the trellis are used.
    pub fn reverse_trellis_full(e_max: usize, n_max: usize, ask: &ASK) -> Vec<Vec<Integer>> {
        let max_energy_step = (e_max - n_max) / 8;
        let mut reverse_trellis = vec![vec![Integer::from(0); max_energy_step + 1]; n_max + 1];
        reverse_trellis[0][0] = Integer::from(1);

        for n in 0..n_max {
            for e in Energies::for_trellis_level(n, n_max, e_max) {
                for amplitude in Amplitudes::possible_from_trellis_node(n, e, n_max, e_max, ask) {
                    let e_next = e + amplitude * amplitude;
                    let current_value = Trellis::get(&reverse_trellis, n, e);
                    Trellis::add(&mut reverse_trellis, n + 1, e_next, current_value);
                }
            }
        }
        reverse_trellis
    }

    /// Calculate the reverse trellis for a given first lexically abandoned sequence (FLAS)
    ///
    /// A number at location (`n`, `e`) in the reverse trellis is the number of length `n` sequences with
    /// the energy `e`. The effect of not using sequences lexically above (and including) the FLAS
    /// is taken into account. Most importantly a number at (`n_max`, `e`) signifies the number of
    /// sequences with energy `e` in an energy trellis with the same parameters as this function,
    /// if only sequences lexically strictly below the FLAS are used.
    ///
    /// The returned reverse trellis is in the shape of a 2D [Vec]. To access values via their `n`
    /// and `e` indexes use the [Trellis::get()] function.
    ///
    /// The calculation is based on formula (20) in section III-C of
    /// <https://doi.org/10.1109/JLT.2022.3201901>.
    pub fn reverse_trellis(
        e_max: usize,
        n_max: usize,
        ask: &ASK,
        first_abandoned_sequence: &[usize],
    ) -> Vec<Vec<Integer>> {
        let max_energy_step = (e_max - n_max) / 8;
        let mut reverse_trellis = vec![vec![Integer::from(0); max_energy_step + 1]; n_max + 1];

        for n in 0..n_max {
            let abandoned_energy = first_abandoned_sequence[0..n]
                .iter()
                .fold(0, |e_total, amplitude| e_total + amplitude * amplitude);
            for e in Energies::for_trellis_level(n, n_max, e_max) {
                for amplitude in Amplitudes::possible_from_trellis_node(n, e, n_max, e_max, ask) {
                    let e_next = e + amplitude * amplitude;
                    let current_value = Trellis::get(&reverse_trellis, n, e);
                    let lower_sequence =
                        i32::from(e == abandoned_energy && amplitude < first_abandoned_sequence[n]);
                    let value_to_add = current_value + lower_sequence;
                    Trellis::add(&mut reverse_trellis, n + 1, e_next, value_to_add);
                }
            }
        }
        reverse_trellis
    }
    /// Get function for trellis values using `n` and `e` as index
    ///
    /// This method translates an index in the form (`n`, `e`) to the correct index in a
    /// 2D [Vec] used for storing the trellis. It then returns the corresponding
    /// value from the provided 2D [Vec] `data`, which is assumed to store a trellis.
    ///
    /// Panics if an invalid `e` value is provided
    pub fn get(data: &[Vec<Integer>], n: usize, e: usize) -> Integer {
        let e_step = (e - n) / 8;
        assert_eq!(e_step * 8, e - n);
        //println!("get: n = {}, e = {} -> {}", n, e, &data[n][e_step]);
        data[n][e_step].clone()
    }
    /// Set function for trellis values using `n` and `e` as index
    fn set(data: &mut [Vec<Integer>], n: usize, e: usize, value: Integer) {
        //println!("set: n = {}, e = {} ({})", n, e, &value);
        let e_step = (e - n) / 8;
        data[n][e_step] = value;
    }
    /// Function to add a value to a trellis value using `n` and `e` as index
    fn add(data: &mut [Vec<Integer>], n: usize, e: usize, value: Integer) {
        //println!("set: n = {}, e = {} ({})", n, e, &value);
        let e_step = (e - n) / 8;
        data[n][e_step] += value;
    }
}

impl Trellis {
    /// Getter method for the values in this trellis
    pub fn value(&self, n: usize, e: usize) -> Integer {
        Trellis::get(&self.data, n, e)
    }
    /// Panics if the amplitude sequence is invalid.
    pub fn validate_amplitudes(&self, amplitude_sequence: &[usize]) -> Result<(), &str> {
        if amplitude_sequence.len() != self.n_max {
            Err("Lenght of amplitude sequence does not equal the length for which the trellis was initialized!")
        } else if amplitude_sequence
            .iter()
            .fold(0, |total, a| cmp::max(total, *a))
            > self.ask.max_amplitude()
        {
            Err("The amplitude sequence contains amplitudes not in the ASK that was used to create the trellis!")
        } else if amplitude_sequence.iter().any(|a| (a % 2 == 0)) {
            Err("All amplitude values must be odd numbers!")
        } else if amplitude_sequence.iter().fold(0, |total, a| total + a * a) > self.e_max {
            Err("The total energy of the sequence exceeds the maximum energy used to calculate the trellis!")
        } else {
            Ok(())
        }
    }
    /// Panics if the sequence index is invalid
    pub fn validate_index(&self, index: &Integer) -> Result<(), &str> {
        if index >= &self.value(0, 0) {
            Err("The given sequence index is out of bounds!")
        } else {
            Ok(())
        }
    }
    /// Returns the probabilities of the amplitudes as a vector
    ///
    /// The first probability in the vector corresponds to the lowest amplitude.
    ///
    /// WARNING: Effects by limiting the used indexes to a power of two are not regarded!!!
    pub fn amplitude_distribution(&self) -> Vec<f32> {
        let mut amplitude_distribution: Vec<f32> =
            Amplitudes::possible_from_trellis_node(0, 0, self.n_max, self.e_max, &self.ask)
                .map(|a| Rational::from((self.value(1, a * a), self.value(0, 0))).to_f32())
                .collect();

        // append zeros to mach the length of the output to the number of amplitudes
        let len_amplitude_distribution = amplitude_distribution.len();
        let num_amplitudes = self.ask.number / 2;
        if len_amplitude_distribution != num_amplitudes {
            amplitude_distribution
                .append(&mut vec![0.0; num_amplitudes - len_amplitude_distribution]);
        }
        amplitude_distribution
    }
    /// Returns the total number of possible sequences
    pub fn num_sequences(&self) -> Integer {
        self.value(0, 0)
    }
    /// Returns the average energy of all possible amplitude sequences
    ///
    /// WARNING: Effects by limiting the used indexes to a power of two are not regarded!!!
    pub fn average_energy(&self) -> f32 {
        let mut average_symbol_energy = 0.0;
        for (probability, amplitude) in self
            .amplitude_distribution()
            .iter()
            .zip(Amplitudes::all_in(&self.ask))
        {
            average_symbol_energy += probability * (amplitude * amplitude) as f32;
        }
        average_symbol_energy * self.n_max as f32
    }

    /// A reverse trellis for this trellis considering only indexes below the largest power of 2 are used
    ///
    /// For details on the returned reverse trellis see the documentation of [Trellis::reverse_trellis()].
    pub fn reverse_trellis_pow_2(&self) -> Vec<Vec<Integer>> {
        let num_bits = self.num_sequences().significant_bits() - 1;
        let mut power_of_two_index = Integer::from(0);
        power_of_two_index.toggle_bit(num_bits);

        if power_of_two_index == self.num_sequences() {
            Trellis::reverse_trellis_full(self.e_max, self.n_max, &self.ask)
        } else {
            let first_abandoned_sequence = self
                .sequence_for_index(&power_of_two_index)
                .expect("If condition guarantees index exists");
            Trellis::reverse_trellis(self.e_max, self.n_max, &self.ask, &first_abandoned_sequence)
        }
    }

    /// A reverse trellis for this trellis if all indexes are used
    ///
    /// This is a short hand for calling [Trellis::reverse_trellis_full()] with the same arguments
    /// used to construct this [Trellis].
    ///
    /// For details on the returned reverse trellis see the documentation of [Trellis::reverse_trellis()].
    pub fn reverse_trellis_full_(&self) -> Vec<Vec<Integer>> {
        Trellis::reverse_trellis_full(self.e_max, self.n_max, &self.ask)
    }

    /// Returns the amplitude sequence for a given index
    ///
    /// Calculations based on algorithm 1 in section III-C of <https://doi.org/10.1109/TWC.2019.2951139>.
    pub fn sequence_for_index(&self, sequence_index: &Integer) -> Result<Vec<usize>, &str> {
        self.validate_index(sequence_index)?;

        let mut amplitude_sequence = Vec::new();

        let mut cumulated_energy = 0;
        let mut num_sequences_left_below = Integer::from(0);
        for n in 0..self.n_max {
            for a in Amplitudes::upto_including_energy(self.e_max - cumulated_energy, &self.ask) {
                let amplitude_energy = a * a;

                let num_possible_sequences_with_a =
                    self.value(n + 1, cumulated_energy + amplitude_energy);
                // it is impossible to leave all sequences possible with amplitude a below
                let index_just_unreachable_with_a =
                    num_sequences_left_below.clone() + num_possible_sequences_with_a.clone();

                if sequence_index < &index_just_unreachable_with_a {
                    // we can reach the target index using amplitude a
                    amplitude_sequence.push(a);
                    cumulated_energy += amplitude_energy;
                    break;
                } else {
                    // we have to use a higher amplitude to leave below all sequences
                    // possible with a in order reach the target index
                    num_sequences_left_below += num_possible_sequences_with_a
                }
            }
        }
        Ok(amplitude_sequence)
    }

    /// Returns the index for a given amplitude sequence
    ///
    /// Calculations based on algorithm 2 in section III-C of <https://doi.org/10.1109/TWC.2019.2951139>.
    pub fn index_for_sequence(&self, amplitude_sequence: &[usize]) -> Result<Integer, &str> {
        self.validate_amplitudes(amplitude_sequence)?;

        // the index of the sequence, before the number of lower sequences is added
        let mut sequence_index = Integer::from(0);

        // pre compute the total energy of all amplitudes upto and excluding a given index
        let cumulated_energy = amplitude_sequence.iter().fold(vec![0], |mut acc, a| {
            acc.push(a * *a + acc[acc.len() - 1]);
            acc
        });

        // add number of lower sequences to the index
        for amplitude_index in 0..self.n_max {
            let current_amplitude = amplitude_sequence[amplitude_index];

            // sum number of possible sequences where the current amplitude would
            // be lower than the real current amplitude
            for lower_amplitude in Amplitudes::upto(current_amplitude, &self.ask) {
                let next_energy = lower_amplitude * lower_amplitude;
                let next_cumulated_energy = cumulated_energy[amplitude_index] + next_energy;
                let next_amplitude_index = amplitude_index + 1;

                sequence_index += self.value(next_amplitude_index, next_cumulated_energy);
            }
        }
        Ok(sequence_index)
    }
}
