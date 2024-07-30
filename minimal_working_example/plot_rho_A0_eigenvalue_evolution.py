from read_data_scripts import get_Posix_for_quantities, open_Posix_dict

import matplotlib.pyplot as plt
import numpy as np

fig, (ax_rho, ax_A0_induced, ax_omega) = plt.subplots(3)

axes_dict = {"rho":ax_rho, "A0_induced":ax_A0_induced}

m = 0
a = 1

directory = "ambjorn"
desired_quantities = ["rho", "A0_induced"]
posix_dict = get_Posix_for_quantities(m, a, directory=directory)
solution_family_array = open_Posix_dict(posix_dict, desired_quantities=["eigenvalues"] + desired_quantities)

lambda_value_array = solution_family_array["lambda_value"]

for quantity in desired_quantities:
    axes_dict[quantity].set_xlabel("z")
    for i, lambda_value in enumerate(lambda_value_array):
        recovered_quantity = solution_family_array[quantity][i]

        n_points = len(recovered_quantity)
        z = np.linspace(0, 1, n_points) 

        axes_dict[quantity].plot(
                z,
                recovered_quantity,
                'b',
                alpha=0.2 + 0.8 * ((i+1)/len(lambda_value_array))**3
                )

ax_omega.plot(lambda_value_array, solution_family_array["eigenvalues"], 'b')

ax_rho.set_ylabel(r"$\rho$")
ax_A0_induced.set_ylabel(r"$A_0$")
ax_omega.set_xlabel(r"$\lambda$")
ax_omega.set_ylabel(r"$\omega_n$")

fig.tight_layout()
plt.show()
