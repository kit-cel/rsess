mod iterators;
pub mod trellis;

#[cfg(test)]
mod tests;

use std::cmp::Ordering;

use iterators::{Amplitudes, Energies};
use rug::{Integer, Rational};
use trellis::Trellis;

/// Representation of the amplitude shift keying coding scheme
#[derive(Clone)]
pub struct ASK {
    number: usize, // is number the right name???
}

impl ASK {
    /// Returns a new [ASK]
    pub fn new(number: usize) -> ASK {
        // number must be a power of 2
        assert_eq!(number.count_ones(), 1);
        ASK { number }
    }
    /// Returns the maximum amplitude in the ASK
    pub fn max_amplitude(&self) -> usize {
        self.number - 1
    }
}

/// Trait with user interface methods common to distribution matching algorithms like ESS or OESS
pub trait DistributionMatcher {
    /// Returns the number of amplitudes in an amplitude sequence
    fn n_max(&self) -> usize;

    /// Returns the maximum allowed energy for an amplitude sequence
    fn e_max(&self) -> usize;

    /// Returns the number of the used ASK, e.g. 8 if 8-ASK is used
    fn ask_num(&self) -> usize;

    /// Returns the amplitude sequence for the given data or an error message if the data is imvalid
    fn encode(&self, data: &Integer) -> Result<Vec<usize>, &str>;

    /// Returns the data for the given amplitude sequence or an error message if the data is invalid
    fn decode(&self, amplitude_sequence: &[usize]) -> Result<Integer, &str>;

    /// Returns the maximum number of possible amplitude sequences
    ///
    /// WARNING: Effect of limiting the used indexes to a power of two is not regarded!!!
    /// If it should be, use [DistributionMatcher::num_sequences()] instead.
    fn num_sequences_possible(&self) -> Integer;

    /// Returns the number of used amplitude sequences
    ///
    /// This will always be a power of two to have a fixed number of bits per amplitude sequence.
    fn num_sequences(&self) -> Integer {
        let mut zero = Integer::from(0);
        zero.toggle_bit(self.num_data_bits());
        zero
    }

    /// Returns the number of bits encoded per amplitude sequence
    fn num_data_bits(&self) -> u32 {
        self.num_sequences_possible().significant_bits() - 1
    }

    /// Returns the probabilities of the amplitude values
    ///
    /// The probabilities are returned as an array with the lowest index corresponding to the
    /// lowest amplitude.
    fn amplitude_distribution(&self) -> Vec<f32>;

    /// Returns the  probabilities of the amplitude sequence energies
    ///
    /// The probabilities are returned as an array. The corresponding energy values can be
    /// calculated with `n_max + 8i` (`n_max` is the number of amplitudes per sequence and `i`
    /// is the index in the returned array).
    fn energy_distribution(&self) -> Vec<f32>;

    /// Returns the average energy of amplitude sequences
    fn average_energy(&self) -> f32 {
        let amplitude_distribution = self.amplitude_distribution();
        let mut amplitude = 1;
        let mut average_energy = 0.0;
        for amplitude_probability in amplitude_distribution {
            average_energy += amplitude_probability * (amplitude * amplitude) as f32;
            amplitude += 2;
        }
        average_energy
    }
}

/// Encoder / decoder using the ESS algorithm
///
/// The implementation follows <https://doi.org/10.1109/TWC.2019.2951139>.
/// Most of the algorithm is implemented in [trellis::Trellis] though.
pub struct ESS {
    trellis: Trellis,
}

impl ESS {
    /// Returns a new [ESS] instance
    pub fn new(e_max: usize, n_max: usize, ask: ASK) -> ESS {
        let trellis = Trellis::new(e_max, n_max, ask);
        ESS { trellis }
    }
}

impl DistributionMatcher for ESS {
    fn n_max(&self) -> usize {
        self.trellis.n_max
    }
    fn e_max(&self) -> usize {
        self.trellis.e_max
    }
    fn ask_num(&self) -> usize {
        self.trellis.ask.max_amplitude() + 1
    }
    fn encode(&self, data: &Integer) -> Result<Vec<usize>, &str> {
        self.trellis.sequence_for_index(data)
    }
    fn decode(&self, amplitude_sequence: &[usize]) -> Result<Integer, &str> {
        self.trellis.index_for_sequence(amplitude_sequence)
    }
    fn num_sequences_possible(&self) -> Integer {
        self.trellis.num_sequences()
    }
    /// Returns the probabilities of the amplitude values
    ///
    /// The probabilities are returned as an array with the lowest index corresponding to the
    /// lowest amplitude.
    ///
    /// Calculation is based on the formulas (21)-(23) in <https://doi.org/10.1109/JLT.2022.3201901>.
    /// The formula (21) in the paper appears to have two mistakes, which are corrected here.
    /// The second summand in (21) should be:
    // $\sum\limits_{x<z_j} \mathbb{I}[x=a] T_j^{x^2 + \sum_{i=1}^{j-1} z_i^2}$
    #[doc = include_str!("./lib-ESS-amplitude_distribution-rendered-formula.doc.svg")]
    fn amplitude_distribution(&self) -> Vec<f32> {
        if self.num_sequences() == self.num_sequences_possible() {
            return self.trellis.amplitude_distribution();
        }
        let n_max = self.trellis.n_max;
        let e_max = self.trellis.e_max;
        let reverse_trellis = self.trellis.reverse_trellis_pow_2();
        let power_of_two_sequence = self
            .encode(&self.num_sequences())
            .expect("the case num_sequences() == num_sequences_possible() was already handled");
        let power_of_two_sequence_cumulated_energy =
            power_of_two_sequence.iter().fold(vec![0], |mut total, a| {
                total.push(total[total.len() - 1] + a * a);
                total
            });

        let amplitude_frequencies: Vec<f32> = Amplitudes::all_in(&self.trellis.ask)
            .map(|amplitude| {
                let aa = amplitude * amplitude;
                let val = (0..n_max).map(|j| {
                    let mut acc = Integer::from(0);
                    if e_max < aa {
                        return 0.0;
                    }
                    for e in Energies::for_trellis_level(j, n_max, e_max) {
                        if e + aa > e_max - n_max + j + 1 {
                            continue;
                        }

                        acc += self.trellis.value(j + 1, e + aa)
                            * Trellis::get(&reverse_trellis, j, e);
                    }
                    if amplitude < power_of_two_sequence[j] {
                        acc += self
                            .trellis
                            .value(j + 1, aa + power_of_two_sequence_cumulated_energy[j]);
                    }
                    if amplitude == power_of_two_sequence[j] {
                        let mut val = Integer::from(0);
                        for n in 0..j + 1 {
                            for x in Amplitudes::upto(power_of_two_sequence[n], &self.trellis.ask) {
                                val += self.trellis.value(
                                    n + 1,
                                    x * x + power_of_two_sequence_cumulated_energy[n],
                                );
                            }
                        }
                        acc = acc + self.num_sequences() - val;
                    }
                    Rational::from((acc, self.num_sequences())).to_f32()
                });

                let sum_amplitude_frequencies_per_location: f32 = val.sum();

                // calculate amplitude frequency and return it
                sum_amplitude_frequencies_per_location / n_max as f32
            })
            .collect();

        amplitude_frequencies
    }
    fn energy_distribution(&self) -> Vec<f32> {
        let reverse_trellis = self.trellis.reverse_trellis_pow_2();
        reverse_trellis[reverse_trellis.len() - 1]
            .iter()
            .map(|count| Rational::from((count, self.num_sequences())).to_f32())
            .collect()
    }
}

/// Encoder / decoder using the OESS algorithm
///
/// The implementation follows <https://doi.org/10.1109/JLT.2022.3201901>.
/// Most of the algorithm is implemented in [trellis::Trellis] though.
pub struct OESS {
    full_trellis: Trellis,
    partial_trellis: Trellis,
    num_low_energy_sequences: Integer,
}

impl OESS {
    /// Returns a new OESS instance
    ///
    /// OESS is sensitive to the maximum energy parameter. If there is a lower maximum energy
    /// that would allow the same number of bits in the index, all sequences with energy higher
    /// than this lower maximum energy will not be used. Thus the OESS algorithm is only defined
    /// for the lowest maximum energy that allows for a certain number of bits in the index. More
    /// specifically the following condition must hold:
    ///
    /// `|S(e_max - 8)| < 2^k ≤ |S(e_max)|`
    ///
    /// - `|S(E)|` is the cardinality of the set of all amplitude sequences with energy less than or equal to `E`
    /// - `k` is a positive integer (equal to the number of bits in the index)
    ///
    /// To find an acceptable `e_max` argument use [OESS::optimal_e_max()].
    pub fn new(e_max: usize, n_max: usize, ask: ASK) -> OESS {
        let full_trellis = Trellis::new(e_max - 8, n_max, ask.clone());
        let partial_trellis = Trellis::new_partial(e_max, n_max, ask);

        // assert the constraint |S(e_max - 8)| < 2^k ≤ |S(e_max)|
        let num_low_energy_sequences = full_trellis.num_sequences();
        let total_num_sequences = &num_low_energy_sequences + partial_trellis.num_sequences();
        assert!(OESS::a_lt_pow2_leq_b(&num_low_energy_sequences, &total_num_sequences),
                "OESS with the given maximum energy does not make sense, try decreasing it. For more info see the documentation of this function.");

        OESS {
            full_trellis,
            partial_trellis,
            num_low_energy_sequences,
        }
    }

    /// Returns the next lower optimal `e_max` value
    ///
    /// An ESS encoder with the provided arguments can encode a certain number of bits. This
    /// function returns the lowest `e_max` parameter so that an ESS encoder with the new
    /// `e_max` still encodes the same number of bits.
    pub fn optimal_e_max(e_max: usize, n_max: usize, ask: &ASK) -> usize {
        let reverse_trellis = Trellis::reverse_trellis_full(e_max, n_max, ask);

        let mut nums_sequences_possible_with_energy_step =
            reverse_trellis[n_max]
                .iter()
                .fold(vec![Integer::from(0)], |mut total, step| {
                    let num_sequences: Integer = (&total[total.len() - 1] + step).into();
                    total.push(num_sequences);
                    total
                });
        // first element only needed for initialization
        nums_sequences_possible_with_energy_step.remove(0);

        // find the energy step so that $|S(E_max - 8)| < 2^k ≤ |S(E_max)|$
        for (energy_step, num_sequences) in nums_sequences_possible_with_energy_step
            .iter()
            .enumerate()
            .rev()
        {
            if OESS::a_lt_pow2_leq_b(
                &nums_sequences_possible_with_energy_step[energy_step - 1],
                num_sequences,
            ) {
                return 8 * energy_step + n_max;
            }
        }
        panic!("OESS::optimal_e_max: This code should never be reached!");
    }

    /// Returns true, if a < 2^k ≤ b
    fn a_lt_pow2_leq_b(a: &Integer, b: &Integer) -> bool {
        let power_of_two: Integer = if a.is_power_of_two() {
            (a + Integer::from(1)).next_power_of_two_ref().into()
        } else {
            a.next_power_of_two_ref().into()
        };
        b >= &power_of_two
    }
}

impl DistributionMatcher for OESS {
    fn n_max(&self) -> usize {
        self.full_trellis.n_max
    }
    fn e_max(&self) -> usize {
        self.partial_trellis.e_max
    }
    fn ask_num(&self) -> usize {
        self.full_trellis.ask.max_amplitude() + 1
    }
    fn encode(&self, data: &Integer) -> Result<Vec<usize>, &str> {
        if data < &self.num_low_energy_sequences {
            self.full_trellis.sequence_for_index(data)
        } else {
            let partial_trellis_index: Integer = (data - &self.num_low_energy_sequences).into();
            self.partial_trellis
                .sequence_for_index(&partial_trellis_index)
        }
    }
    fn decode(&self, amplitude_sequence: &[usize]) -> Result<Integer, &str> {
        let sequence_energy = amplitude_sequence
            .iter()
            .fold(0, |total_energy, amplitude| {
                total_energy + amplitude * amplitude
            });
        if sequence_energy <= self.full_trellis.e_max {
            self.full_trellis.index_for_sequence(amplitude_sequence)
        } else {
            let num_max_energy_sequences = self
                .partial_trellis
                .index_for_sequence(amplitude_sequence)?;
            Ok((&num_max_energy_sequences + &self.num_low_energy_sequences).into())
        }
    }
    fn num_sequences_possible(&self) -> Integer {
        &self.num_low_energy_sequences + self.partial_trellis.num_sequences()
    }
    /// Returns the probabilities of the amplitude values
    ///
    /// The probabilities are returned as an array with the lowest index corresponding to the
    /// lowest amplitude.
    ///
    /// Calculation for the partial trellis is based on formula (36) in section IV-D in
    /// <https://doi.org/10.1109/JLT.2022.3201901>.
    fn amplitude_distribution(&self) -> Vec<f32> {
        if self.num_sequences() == self.num_sequences_possible() {
            let equivalent_trellis =
                Trellis::new(self.e_max(), self.n_max(), ASK::new(self.ask_num()));
            return equivalent_trellis.amplitude_distribution();
        }

        // Calculation of amplitude counts of full trellis

        let num_amplitudes_full =
            self.full_trellis.num_sequences() * Integer::from(self.full_trellis.n_max);
        let amplitude_counts_full: Vec<Integer> = self
            .full_trellis
            .amplitude_distribution()
            .into_iter()
            .map(|p| {
                Integer::from(
                    (Rational::from_f32(p).unwrap() * Rational::from((&num_amplitudes_full, 1)))
                        .round_ref(),
                )
            })
            .collect();

        // Calculation of amplitude counts of partial trellis

        let n_max = self.partial_trellis.n_max;
        let e_max = self.partial_trellis.e_max;
        let ask = &self.partial_trellis.ask;
        let power_of_two_sequence = self
            .encode(&self.num_sequences())
            .expect("the case num_sequences() == num_sequences_possible() was already handled");
        let power_of_two_sequence_cumulated_energy =
            power_of_two_sequence.iter().fold(vec![0], |mut total, a| {
                total.push(total[total.len() - 1] + a * a);
                total
            });

        let amplitude_counts_partial: Vec<Integer> = Amplitudes::all_in(&self.partial_trellis.ask)
            .map(|amplitude| {
                let aa = amplitude * amplitude;
                let amplitude_counts_per_location = (0..n_max).map(|j| {
                    let first_summand = (0..j)
                        .map(|n| {
                            Amplitudes::upto(power_of_two_sequence[n], ask)
                                .map(|x| {
                                    let energy =
                                        aa + x * x + power_of_two_sequence_cumulated_energy[n];

                                    if energy <= e_max - n_max + n + 2 {
                                        self.partial_trellis.value(n + 2, energy)
                                    } else {
                                        Integer::from(0)
                                    }
                                })
                                .sum::<Integer>()
                        })
                        .sum::<Integer>();

                    let second_summand = match amplitude.cmp(&power_of_two_sequence[j]) {
                        Ordering::Greater => Integer::from(0),
                        Ordering::Less => self
                            .partial_trellis
                            .value(j + 1, aa + power_of_two_sequence_cumulated_energy[j]),
                        Ordering::Equal => (j + 1..n_max)
                            .map(|n| {
                                Amplitudes::upto(power_of_two_sequence[n], ask)
                                    .map(|x| {
                                        self.partial_trellis.value(
                                            n + 1,
                                            x * x + power_of_two_sequence_cumulated_energy[n],
                                        )
                                    })
                                    .sum::<Integer>()
                            })
                            .sum::<Integer>(),
                    };

                    first_summand + second_summand
                });

                let amplitude_count: Integer = amplitude_counts_per_location.sum();
                amplitude_count
            })
            .collect();

        // Total amplitude count is the sum of full and partial trellis counts

        let amplitude_counts =
            std::iter::zip(amplitude_counts_full, amplitude_counts_partial).map(|(a, b)| a + b);

        // Convert counts to frequencies by dividing by the number of amplitudes

        let num_amplitudes = self.num_sequences() * n_max;

        // Calculate amplitude frequencies and return them
        amplitude_counts
            .map(|count| Rational::from((count, &num_amplitudes)).to_f32())
            .collect()
    }
    fn energy_distribution(&self) -> Vec<f32> {
        let reverse_trellis_full = self.full_trellis.reverse_trellis_full_();
        // Energy distribution for full trellis is number of sequences with a given energy divided
        // by the total number of sequences.
        let mut energy_distribution = reverse_trellis_full[self.full_trellis.n_max]
            .iter()
            .map(|count| Rational::from((count, self.num_sequences())).to_f32())
            .collect::<Vec<f32>>();

        // Energy distribution for the partial trellis is only one value because all sequences have
        // the maximum energy.
        let count_partial = self.num_sequences() - self.full_trellis.num_sequences();
        let fraction_partial = Rational::from((count_partial, self.num_sequences())).to_f32();

        // Combine results from full and partial trellis
        energy_distribution.push(fraction_partial);

        energy_distribution
    }
}
