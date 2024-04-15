import numpy as np
import matplotlib.pyplot as plt

from scripts.read_files import read_files

m=3.

data = read_files(m=m)

omegas = np.array([np.genfromtxt(filename, delimiter="\n") for filename in data["eigenvalue_array"]])
lambda_value_array = data["lambda_value"]

#for lambda_value, fixed_lambda_omegas in zip(lambda_value_array, omegas):
#    plt.plot(fixed_lambda_omegas, [lambda_value]*len(fixed_lambda_omegas), 'x')

fig, ax = plt.subplots()
for i, omega in enumerate(omegas[0]):
    ax.plot(omegas[:,i], lambda_value_array, 'b')

ax.set_title(f'Evolution of the eigenvalues with $\lambda$, m={m}')
ax.set_xlabel('$\omega_n$')
ax.set_ylabel('$\lambda$')

plt.show()
