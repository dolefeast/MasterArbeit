from read_data_scripts import get_Posix_for_quantities, open_Posix_dict
from minimal_working_example import root_mean_square
from read_files import get_directory_m_a, str_to_float

import matplotlib.pyplot as plt


bcs = "dirichlet"
max_lambda_density=50/16

directory, m, a = get_directory_m_a(bcs)

desired_quantities = ["rho", "A0_induced"]
posix_dict = get_Posix_for_quantities(m, a, directory=directory, max_lambda_density=max_lambda_density)

solution_family_array = open_Posix_dict(posix_dict, desired_quantities=desired_quantities)

lambda_value_array = solution_family_array["lambda_value"]

L2_norms = {
        "rho": [], 
        "A0_induced": [], 
        }

fig = plt.figure(f'{directory}, m={str_to_float(m)}, a={str_to_float(a)}')
(ax_rho, ax_A0_induced) = fig.subplots(2)
axes_dict = {"rho":ax_rho, "A0_induced":ax_A0_induced}

for quantity in desired_quantities:
    axes_dict[quantity].set_xlabel("$\lambda$")
    for i, lambda_value in enumerate(lambda_value_array):
            L2_norms[quantity].append(
                    root_mean_square(
                    solution_family_array[quantity][i]
                    )
                    )
    axes_dict[quantity].plot(lambda_value_array, L2_norms[quantity])

ax_rho.set_ylabel(r"$\sqrt{\langle \rho^2\rangle}$")
ax_A0_induced.set_ylabel(r"$\sqrt{\langle A_0^2\rangle}$")

fig.tight_layout()
plt.show()
