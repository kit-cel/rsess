import pyrsess
import matplotlib.pyplot as plt
import numpy as np

max_energy = 1120
seq_len = 96
max_amplitude = 8
myshaper = pyrsess.ESS(max_energy, seq_len, max_amplitude)

# encode a random bit array
index_len = myshaper.num_data_bits()
print(f'Number of input bits: {myshaper.num_data_bits()}')

rand_bits = np.random.randint(2, size=index_len)
print(f'Random bits to be encoded: {rand_bits}')
print(f'Amplitude sequence: {myshaper.encode(rand_bits)}')

# test encoding and decoding random bit arrays
num_transmissions = 1000
idxs_correct = np.random.randint(0, 2, (num_transmissions, index_len))
results = myshaper.multi_encode(idxs_correct)

idxs = myshaper.multi_decode(results)

print(f'Correct encoding and decoding of {num_transmissions} indexes: ', end='')
print(np.all([idx_correct == idx for idx_correct, idx in zip(idxs_correct, idxs)]))

# print the amplitude distribution and average energy
print('Amplitude distribution:')
print(myshaper.amplitude_distribution())

print('Average energy:')
print(myshaper.average_energy())

# plot the calculated amplitude distribution vs. the empirical one
flat_results = np.array(results).flatten()
amplitudes, counts = np.unique(flat_results, return_counts=True)

fig, ax = plt.subplots(1,1)
ax.set_xticks(amplitudes)
ax.bar(amplitudes-0.4, counts/flat_results.size, color='b', label='Empirical')
ax.bar(amplitudes+0.4, myshaper.amplitude_distribution(), color='r', label='Calculated')
ax.legend()
ax.set_title('Amplitude Distribution: Empirical vs. Calculated')
ax.set_ylabel('Probability')
ax.set_xlabel('Amplitude')

plt.show()
