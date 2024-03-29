import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path

eigenstate_array = list(
    Path("saved_solutions/dirichlet/normalized_eigenstate").glob("*lambda_20_462_mass_3_0*.txt")
)

eigenstate2 = np.genfromtxt(eigenstate_array[0], delimiter=",")
#eigenstate3 = np.genfromtxt(eigenstate_array[3], delimiter=",")

N_modes = len(eigenstate2)
for i, eigenstate_file in enumerate(eigenstate_array):
    alpha = 0.2 + 0.8 * ((i + 1) / len(eigenstate_array)) ** 3
    #eigenstate = np.genfromtxt(eigenstate_file, delimiter=",")[N_modes // 2 + 3]
    eigenstate = np.genfromtxt(eigenstate_file, delimiter=",")[0]
    print(eigenstate[1])
    plt.plot(
        eigenstate,
        'kx',
        alpha=alpha
    )


plt.show()
