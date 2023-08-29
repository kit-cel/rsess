use rug::Integer;
use rug::Rational;

use crate::trellis::Trellis;
use crate::{DistributionMatcher, ASK, ESS, OESS};

// Trellis
// ------------------------------------------------------

#[test]
fn trellis_paper_example() {
    let trellis = Trellis::new(28, 4, ASK::new(8));
    assert_eq!(trellis.value(0, 0), 19);
    assert_eq!(trellis.value(1, 1), 11);
    assert_eq!(trellis.value(1, 9), 7);
    assert_eq!(trellis.value(1, 25), 1);
    assert_eq!(trellis.value(2, 2), 6);
    assert_eq!(trellis.value(2, 10), 4);
    assert_eq!(trellis.value(2, 18), 3);
    assert_eq!(trellis.value(2, 26), 1);
    assert_eq!(trellis.value(3, 3), 3);
    assert_eq!(trellis.value(3, 11), 2);
    assert_eq!(trellis.value(3, 19), 2);
    assert_eq!(trellis.value(3, 27), 1);
    assert_eq!(trellis.value(4, 4), 1);
    assert_eq!(trellis.value(4, 12), 1);
    assert_eq!(trellis.value(4, 20), 1);
    assert_eq!(trellis.value(4, 28), 1);
}

#[test]
fn trellis_amplitude_distribution() {
    let trellis = Trellis::new(120, 8, ASK::new(8));

    let mut amplitude_distribution_direct = vec![0, 0, 0, 0];

    for index in 0..trellis.num_sequences().to_usize().unwrap_or(0) {
        for amplitude in trellis
            .sequence_for_index(&Integer::from(index))
            .expect("static")
        {
            amplitude_distribution_direct[(amplitude / 2) as usize] += 1;
        }
    }

    let total_amplitudes = trellis.num_sequences().to_f32() * trellis.n_max as f32;
    let amplitude_distribution_trellis: Vec<u32> = trellis
        .amplitude_distribution()
        .iter()
        .map(|p| (total_amplitudes * p).round() as u32)
        .collect();

    assert_eq!(
        amplitude_distribution_direct,
        amplitude_distribution_trellis
    );
}

#[test]
fn trellis_average_energy() {
    let trellis = Trellis::new(80, 5, ASK::new(8));

    let mut average_energy = Rational::from(0);
    for index in 0..trellis.num_sequences().to_u32().unwrap_or(0) {
        let energy = trellis
            .sequence_for_index(&Integer::from(index))
            .expect("static")
            .iter()
            .fold(0, |total, a| total + a * a);
        average_energy += Rational::from((energy, trellis.num_sequences()));
    }
    let average_energy = average_energy.to_f32();

    assert_eq!(average_energy, trellis.average_energy());
}

#[test]
fn trellis_reverse_trellis_paper_example() {
    let reverse_trellis = Trellis::reverse_trellis(60, 4, &ASK::new(8), &vec![5, 1, 3, 1]);
    let reverse_trellis_from_paper = vec![
        vec![0, 0, 0, 0, 0, 0, 0, 0],
        vec![1, 1, 0, 0, 0, 0, 0, 0],
        vec![1, 2, 1, 1, 1, 0, 1, 1],
        vec![1, 3, 3, 4, 4, 2, 3, 5],
        vec![1, 4, 6, 8, 11, 9, 10, 15],
    ];
    assert_eq!(reverse_trellis, reverse_trellis_from_paper);
}

#[test]
fn trellis_reverse_trellis_full() {
    let e_max = 2500;
    let n_max = 128;
    let reverse_trellis = Trellis::reverse_trellis_full(e_max, n_max, &ASK::new(8));
    let trellis = Trellis::new(e_max, n_max, ASK::new(8));
    assert_eq!(
        trellis.num_sequences(),
        reverse_trellis[n_max].iter().sum::<Integer>()
    );

    let e_max = 7 * 7 * 96;
    let n_max = 96;
    let reverse_trellis = Trellis::reverse_trellis_full(e_max, n_max, &ASK::new(8));
    let trellis = Trellis::new(e_max, n_max, ASK::new(8));
    assert_eq!(
        trellis.num_sequences(),
        reverse_trellis[n_max].iter().sum::<Integer>()
    );
}

// ESS
// ------------------------------------------------------

#[test]
fn ess_static_shape_deshape_all_sequences() {
    // trellis with effectively no energy constraint
    let ess = ESS::new(200, 3, ASK::new(8));
    assert_eq!(ess.num_sequences_possible(), 64);
    for index in 0..64 {
        let amplitudes = [1, 3, 5, 7];
        // because of no energy constraints, constructing the sequence
        // becomes converting binary to base 4
        let sequence = vec![
            amplitudes[(0b110000 & index) >> 4],
            amplitudes[(0b001100 & index) >> 2],
            amplitudes[0b000011 & index],
        ];
        assert_eq!(ess.decode(&sequence).expect("static"), index as u32);
        assert_eq!(ess.encode(&Integer::from(index)).expect("static"), sequence);
    }
}

#[test]
fn ess_dynamic_shape_deshape() {
    let ess = ESS::new(1120, 96, ASK::new(8));
    let num_sequences = ess.num_sequences_possible();

    let mut index = Integer::from(1);

    while index <= num_sequences {
        let sequence = ess.encode(&index).expect("static");
        assert_eq!(ess.decode(&sequence).expect("static"), index);

        index = &index * Integer::from(2);
    }
}

#[test]
fn ess_dynamic_shape_deshape_large() {
    let ess = ESS::new(20_000, 200, ASK::new(16));
    let num_sequences = ess.num_sequences_possible();

    let mut index = Integer::from(1);

    while index <= num_sequences {
        let sequence = ess.encode(&index).expect("static");
        assert_eq!(ess.decode(&sequence).expect("static"), index);

        index = &index * Integer::from(20);
    }
}

#[test]
fn ess_amplitude_distribution() {
    fn test_distribution(e_max: usize, n_max: usize, ask: ASK) {
        println!(
            "test_distribution({}, {}, ASK::new({}))",
            e_max, n_max, ask.number
        );
        let ess = ESS::new(e_max, n_max, ask);
        let mut distribution = vec![0, 0, 0, 0];
        for index in 0..ess.num_sequences().to_u32().unwrap_or(0) {
            let sequence = ess.encode(&index.into()).expect("static");
            for amplitude in sequence {
                match amplitude {
                    1 => distribution[0] += 1,
                    3 => distribution[1] += 1,
                    5 => distribution[2] += 1,
                    7 => distribution[3] += 1,
                    _ => (),
                }
            }
        }
        let num_amplitudes = n_max as f32 * ess.num_sequences().to_f32();
        let ess_distribution: Vec<u32> = ess
            .amplitude_distribution()
            .into_iter()
            .map(|val| (num_amplitudes * val).round() as u32)
            .collect();
        assert_eq!(ess_distribution, distribution);
    }

    test_distribution(28, 4, ASK::new(8));
    test_distribution(60, 4, ASK::new(8));
    test_distribution(7 * 7 * 4, 4, ASK::new(8));
    test_distribution(196, 4, ASK::new(8));
    test_distribution(200, 8, ASK::new(8));
    test_distribution(100, 6, ASK::new(8));
}

#[test]
fn ess_energy_distribution() {
    let e_max = 60;
    let n_max = 4;
    let ask = ASK::new(8);
    let ess = ESS::new(e_max, n_max, ask);

    let mut distribution = vec![0; (e_max - n_max) / 8 + 1];
    for index in 0..ess.num_sequences().to_u32().unwrap_or(0) {
        let sequence = ess.encode(&index.into()).expect("static");
        let sequence_energy = sequence.into_iter().fold(0, |total, new| total + new * new);
        let l = (sequence_energy - n_max) / 8;
        distribution[l] += 1;
    }

    let ess_distribution: Vec<u32> = ess
        .energy_distribution()
        .into_iter()
        .map(|val| (ess.num_sequences().to_f32() * val).round() as u32)
        .collect();

    assert_eq!(ess_distribution, distribution);
}

// OESS
// ------------------------------------------------------

#[test]
fn oess_static_shape_deshape_all_sequences() {
    let oess = OESS::new(28, 4, ASK::new(8));
    let sequences = [
        vec![1, 1, 1, 1],
        vec![1, 1, 1, 3],
        vec![1, 1, 3, 1],
        vec![1, 1, 3, 3],
        vec![1, 3, 1, 1],
        vec![1, 3, 1, 3],
        vec![1, 3, 3, 1],
        vec![3, 1, 1, 1],
        vec![3, 1, 1, 3],
        vec![3, 1, 3, 1],
        vec![3, 3, 1, 1],
        vec![1, 1, 1, 5],
        vec![1, 1, 5, 1],
        vec![1, 3, 3, 3],
        vec![1, 5, 1, 1],
        vec![3, 1, 3, 3],
        vec![3, 3, 1, 3],
        vec![3, 3, 3, 1],
        vec![5, 1, 1, 1],
    ];

    assert_eq!(oess.full_trellis.num_sequences(), 11);
    assert_eq!(oess.partial_trellis.num_sequences(), 8);

    for (index, sequence) in sequences.iter().enumerate() {
        assert_eq!(oess.decode(&sequence).expect("static"), index);
        assert_eq!(
            &oess.encode(&Integer::from(index)).expect("static"),
            sequence
        );
    }
}

#[test]
fn oess_lexicographic_ordering() {
    let ask_num = 8;
    let e_max = 60;
    let oess = OESS::new(e_max, 4, ASK::new(ask_num));
    let mut last_order_num = 0;
    let mut last_order_num_e_max = 0;
    let mut had_e_max = false;

    for index in 0..(oess.num_sequences().to_u32().unwrap_or(0)) {
        let sequence = oess.encode(&Integer::from(index)).expect("static");

        let energy: usize = (&sequence).into_iter().map(|a| a * a).sum();
        let order_num = (&sequence)
            .into_iter()
            .rev()
            .enumerate()
            .fold(0, |total, (idx, a)| total + a * ask_num.pow(idx as u32));

        println!("sequence: {:?}, energy: {}", sequence, energy);
        println!(
            "order_num: {}, last_order_num: {}",
            order_num, last_order_num
        );
        if energy == e_max {
            assert!(order_num > last_order_num_e_max);
            last_order_num_e_max = order_num;
            had_e_max = true;
        } else {
            assert!(order_num > last_order_num && !had_e_max);
            last_order_num = order_num;
        }
    }
}

#[test]
fn oess_dynamic_shape_deshape() {
    let oess = OESS::new(1120, 96, ASK::new(8));
    let num_sequences = oess.num_sequences_possible();

    let mut index = Integer::from(1);

    while index <= num_sequences {
        let sequence = oess.encode(&index).expect("static");
        assert_eq!(oess.decode(&sequence).expect("static"), index);

        index = &index * Integer::from(2);
    }
}

#[test]
fn oess_dynamic_shape_deshape_large() {
    let oess = OESS::new(16_992, 200, ASK::new(16));
    let num_sequences = oess.num_sequences_possible();

    let mut index = Integer::from(1);

    while index <= num_sequences {
        let sequence = oess.encode(&index).expect("static");
        assert_eq!(oess.decode(&sequence).expect("static"), index);

        index = &index * Integer::from(20);
    }
}

#[test]
fn oess_find_optimal_energy() {
    let e_max = 120;
    let n_max = 8;
    let ask = ASK::new(8);
    let optimal_e_max = OESS::optimal_e_max(e_max, n_max, &ask);
    OESS::new(optimal_e_max, n_max, ask); // if energy is optimal this does not panic

    let e_max = 1008;
    let n_max = 128;
    let ask = ASK::new(8);
    let optimal_e_max = OESS::optimal_e_max(e_max, n_max, &ask);
    OESS::new(optimal_e_max, n_max, ask); // if energy is optimal this does not panic

    let n_max = 128;
    let ask = ASK::new(8);
    let e_max = 7 * 7 * n_max;
    let optimal_e_max = OESS::optimal_e_max(e_max, n_max, &ask);
    assert_eq!(optimal_e_max, e_max);
    OESS::new(optimal_e_max, n_max, ask); // if energy is optimal this does not panic
}

#[test]
fn oess_amplitude_distribution() {
    fn test_distribution(e_max: usize, n_max: usize, ask: ASK) {
        println!(
            "test_distribution({}, {}, ASK::new({}))",
            e_max, n_max, ask.number
        );
        let oess = OESS::new(e_max, n_max, ask);
        let mut distribution = vec![0, 0, 0, 0];
        for index in 0..oess.num_sequences().to_u32().unwrap_or(0) {
            let sequence = oess.encode(&index.into()).expect("static");
            for amplitude in sequence {
                match amplitude {
                    1 => distribution[0] += 1,
                    3 => distribution[1] += 1,
                    5 => distribution[2] += 1,
                    7 => distribution[3] += 1,
                    _ => (),
                }
            }
        }
        let num_amplitudes = n_max as f32 * oess.num_sequences().to_f32();
        let oess_distribution: Vec<u32> = oess
            .amplitude_distribution()
            .into_iter()
            .map(|val| (num_amplitudes * val).round() as u32)
            .collect();
        assert_eq!(oess_distribution, distribution);
    }

    test_distribution(28, 4, ASK::new(8));
    test_distribution(60, 4, ASK::new(8));
    test_distribution(7 * 7 * 4, 4, ASK::new(8));
    test_distribution(112, 8, ASK::new(8));
    test_distribution(100, 6, ASK::new(8));
}

#[test]
fn oess_energy_distribution() {
    let e_max = 60;
    let n_max = 4;
    let ask = ASK::new(8);
    let oess = OESS::new(e_max, n_max, ask);

    let mut distribution = vec![0; (e_max - n_max) / 8 + 1];
    for index in 0..oess.num_sequences().to_u32().unwrap_or(0) {
        let sequence = oess.encode(&index.into()).expect("static");
        let sequence_energy = sequence.into_iter().fold(0, |total, new| total + new * new);
        let l = (sequence_energy - n_max) / 8;
        distribution[l] += 1;
    }

    let oess_distribution: Vec<u32> = oess
        .energy_distribution()
        .into_iter()
        .map(|val| (oess.num_sequences().to_f32() * val).round() as u32)
        .collect();

    assert_eq!(oess_distribution, distribution);
}
