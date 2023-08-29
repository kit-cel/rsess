use crate::ASK;

/// An iterator to iterate over amplitudes
pub struct Amplitudes {
    current: usize,
    max: usize,
    including: bool,
}

impl Amplitudes {
    /// Returns an iterator over amplitudes up to and including a given max energy
    ///
    /// The first amplitude is 1 and the last amplitude is the largest amplitude with energy
    /// less than or equal to `max_energy`. The Iterator only returns amplitudes available in the
    /// given ASK.
    ///
    /// Examples (assuming 8-ASK):
    /// - `max_energy = 25`: Iterator amplitudes: 1, 3, 5
    /// - `max_energy = 90`: Iterator amplitudes: 1, 3, 5, 7 (9 is not in the iterator because it
    /// is not in an 8-ASK)
    pub fn upto_including_energy(max_energy: usize, ask: &ASK) -> Amplitudes {
        let max_amplitude = f64::floor(f64::sqrt(max_energy as f64)) as usize;
        Amplitudes {
            current: 1,
            max: std::cmp::min(max_amplitude, ask.max_amplitude()),
            including: true,
        }
    }
    /// Returns an iterator over all amplitudes that are possible from a given trellis node
    ///
    /// An amplitude is considered impossible, if it makes it impossible for the full amplitude
    /// sequence to have a maximum energy less than or equal to `e_max`.
    pub fn possible_from_trellis_node(
        n: usize,
        e: usize,
        n_max: usize,
        e_max: usize,
        ask: &ASK,
    ) -> Amplitudes {
        // for n+1 as we must consider the energy constraint of the next trellis node
        let next_local_max_energy = e_max - n_max + n + 1;
        let remaining_energy = next_local_max_energy - e;
        let max_amplitude = f64::floor(f64::sqrt(remaining_energy as f64)) as usize;
        Amplitudes {
            current: 1,
            max: std::cmp::min(max_amplitude, ask.max_amplitude()),
            including: true,
        }
    }
    /// Returns an iterator over all amplitudes upto and excluding the given `max_amplitude`
    ///
    /// If `max_amplitude` exeeds the amplitudes available in the used ASK, the iterator ends with
    /// the highest amplitude available in the ASK.
    pub fn upto(max_amplitude: usize, ask: &ASK) -> Amplitudes {
        Amplitudes {
            current: 1,
            max: std::cmp::min(max_amplitude, ask.max_amplitude() + 1),
            including: false,
        }
    }
    /// Returns an iterator over all amplitudes in the given ASK
    pub fn all_in(ask: &ASK) -> Amplitudes {
        Amplitudes {
            current: 1,
            max: ask.max_amplitude(),
            including: true,
        }
    }
}

impl Iterator for Amplitudes {
    type Item = usize;

    /// Returns the next amplitude
    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.max || (self.including && self.current == self.max) {
            let return_amplitude = self.current;
            self.current += 2;

            Some(return_amplitude)
        } else {
            None
        }
    }
}

/// An iterator over energy levels in the energy trellis
pub struct Energies {
    current: usize,
    max: usize,
    including: bool,
}

impl Energies {
    /// Returns an iterator over all energy levels at location `n` in the energy trellis
    pub fn for_trellis_level(n: usize, n_max: usize, e_max: usize) -> Energies {
        let local_max_energy = if e_max + n >= n_max {
            e_max + n - n_max
        } else {
            0
        };
        Energies {
            current: n,
            max: local_max_energy,
            including: true,
        }
    }
}

impl Iterator for Energies {
    type Item = usize;

    /// Returns the next energy level
    fn next(&mut self) -> Option<Self::Item> {
        if self.current < self.max || (self.including && self.current == self.max) {
            let return_energy = self.current;
            self.current += 8;

            Some(return_energy)
        } else {
            None
        }
    }
}
