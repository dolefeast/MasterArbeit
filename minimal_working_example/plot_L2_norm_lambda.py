from read_data_scripts import get_Posix_for_quantities, open_Posix_dict
from minimal_working_example import root_mean_square

import matplotlib.pyplot as plt

fig, (ax_rho, ax_A0_induced) = plt.subplots(2)

axes_dict = {"rho":ax_rho, "A0_induced":ax_A0_induced}

m = 0
a = 1

directory = "ambjorn"
desired_quantities = ["rho", "A0_induced"]
posix_dict = get_Posix_for_quantities(m, a, directory=directory)
solution_family_array = open_Posix_dict(posix_dict, desired_quantities=desired_quantities)

lambda_value_array = solution_family_array["lambda_value"]

for quantity in desired_quantities:
    axes_dict[quantity].set_xlabel("z")
    for i, lambda_value in enumerate(lambda_value_array):
        axes_dict[quantity].plot(
                lambda_value,
                root_mean_square(
                    solution_family_array[quantity][i]
                    ),
                'bo'
                )

ax_rho.set_ylabel(r"$\rho$")
ax_A0_induced.set_ylabel(r"$A_0$")

fig.tight_layout()
plt.show()
