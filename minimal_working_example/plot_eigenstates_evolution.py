from read_data_scripts import get_Posix_for_quantities, open_Posix_dict

import matplotlib.pyplot as plt
import numpy as np


fig, axes_list = plt.subplots(4, 2, figsize=(16, 9))

#len(axes_list) = 4, len(axes_list[0]) = 2

modes = [-2, -1, 0, 1]
m = 0
a = 1

directory = "ambjorn"
desired_quantities = ["eigenvalues", "eigenstates", "A0_induced"]
posix_dict = get_Posix_for_quantities(m, a, directory=directory)
solution_family_array = open_Posix_dict(posix_dict, desired_quantities=desired_quantities)

lambda_value_array = solution_family_array["lambda_value"]

for i, (eigenvalue_array, eigenstate_array, A0_induced, lambda_value) in enumerate(zip(
    solution_family_array["eigenvalues"],
    solution_family_array["eigenstates"],
    solution_family_array["A0_induced"],
    solution_family_array["lambda_value"],
    )):
        for mode, (ax_mode, ax_rho_mode) in zip(modes, axes_list):
            max_N = len(eigenvalue_array) // 2 # The cutoff N
                
            mode_index = max_N+mode
            if mode_index < 0 or mode_index>max_N:
                continue

            eigenstate = eigenstate_array[mode_index]
            eigenvalue = eigenvalue_array[mode_index]

            n_points = len(eigenstate)
            z = np.linspace(0, 1, n_points) 

            ax_mode.plot(
                    z,
                    eigenstate,
                    'b',
                    alpha=0.2 + 0.8 * ((i+1)/len(lambda_value_array))**3
                    )

            rho_mode = 1/2 * (eigenvalue + lambda_value * (z-1/2) - A0_induced) * eigenstate**2

            ax_rho_mode.plot(
                    z,
                    rho_mode,
                    'b',
                    alpha=0.2 + 0.8 * ((i+1)/len(lambda_value_array))**3
                    )

            ax_mode.set_xlabel(r"z")
            ax_rho_mode.set_xlabel(r"z")
            ax_mode.set_ylabel(r"$\phi_{"+ "{}".format(mode + int(mode>=0))+"}$")
            ax_rho_mode.set_ylabel(r"$\rho_{"+ "{}".format(mode + int(mode>=0))+"}$")

            

fig.tight_layout()
plt.show()
