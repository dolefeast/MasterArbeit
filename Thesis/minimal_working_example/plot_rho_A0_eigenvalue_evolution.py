from read_data_scripts import get_Posix_for_quantities, open_Posix_dict
from read_files import get_directory_m_a, str_to_float

import matplotlib.pyplot as plt
import numpy as np


bcs = "dirichlet"
max_lambda_density = None

directory, m, a = get_directory_m_a(bcs)
desired_quantities = ["rho", "A0_induced"]
posix_dict = get_Posix_for_quantities(m, a, directory=directory, max_lambda_density=max_lambda_density)
solution_family_array = open_Posix_dict(posix_dict, desired_quantities=["eigenvalues"] + desired_quantities)

maxN = len(solution_family_array["eigenvalues"][0])//2
if maxN == 1:
    omega_ylim = solution_family_array["eigenvalues"][0][:]
elif maxN >= 2:
    omega_ylim = (-10, 10)

lambda_value_array = solution_family_array["lambda_value"]

print("The maximum lambda value obtained is", lambda_value_array[-1]) #, ", with a minimum lambda step of", lambda_value_array[-1]-lambda_value_array[-2])

fig = plt.figure(f'{directory}, m={str_to_float(m)}, a={str_to_float(a)}')
(ax_rho, ax_A0_induced, ax_omega) = fig.subplots(3)
axes_dict = {"rho":ax_rho, "A0_induced":ax_A0_induced}

for quantity in desired_quantities:
    axes_dict[quantity].set_xlabel("z")
    for i, lambda_value in enumerate(lambda_value_array):
        recovered_quantity = solution_family_array[quantity][i]

        n_points = len(recovered_quantity)
        z = np.linspace(0, 1, n_points) 

        alpha = 0.2 + 0.8 * ((i+1)/len(lambda_value_array))**3

        axes_dict[quantity].plot(
                z,
                recovered_quantity,
                'b',
                alpha=alpha,
                )


ax_omega.plot(lambda_value_array, solution_family_array["eigenvalues"], 'b')

ax_rho.set_ylabel(r"$\rho$")
ax_A0_induced.set_ylabel(r"$A_0$")
ax_omega.set_xlabel(r"$\lambda$")
ax_omega.set_ylabel(r"$\omega_n$")



ax_omega.set_ylim(omega_ylim)

fig.tight_layout()
plt.show()
