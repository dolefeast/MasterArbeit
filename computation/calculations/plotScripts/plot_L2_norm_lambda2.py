from read_data_scripts import get_Posix_for_quantities, open_Posix_dict
from read_files import get_directory_m_a, str_to_float

import matplotlib.pyplot as plt

bcs = "dirichlet"
max_lambda_density=None

directory, m, a = get_directory_m_a(bcs)
directory1, m1, a1 = get_directory_m_a(bcs)
directory2, m1, a1 = get_directory_m_a(bcs)

desired_quantities = ["rho", "A0_induced"]

posix_dict = get_Posix_for_quantities(m, a, directory=directory, max_lambda_density=max_lambda_density)
posix_dict1 = get_Posix_for_quantities(m, a, directory=directory1, max_lambda_density=max_lambda_density)
posix_dict2 = get_Posix_for_quantities(m, a, directory=directory2, max_lambda_density=max_lambda_density)

solution_family_array = open_Posix_dict(posix_dict, desired_quantities=desired_quantities)
solution_family_array1 = open_Posix_dict(posix_dict1, desired_quantities=desired_quantities)
solution_family_array2 = open_Posix_dict(posix_dict2, desired_quantities=desired_quantities)

lambda_value_array = solution_family_array["lambda_value"]
lambda_value_array1 = solution_family_array1["lambda_value"]
lambda_value_array2 = solution_family_array2["lambda_value"]

L2_norms = {
        "rho": [], 
        "A0_induced": [], 
        }

L2_norms1 = {
        "rho": [], 
        "A0_induced": [], 
        }
L2_norms2 = {
        "rho": [], 
        "A0_induced": [], 
        }

fig = plt.figure() #f'{directory}, m={str_to_float(m)}, a={str_to_float(a)}')
(ax_rho, ax_A0_induced) = fig.subplots(2)
axes_dict = {"rho":ax_rho, "A0_induced":ax_A0_induced}

for quantity in desired_quantities:
    axes_dict[quantity].set_xlabel("$\lambda$")
    for i, lambda_value in enumerate(lambda_value_array):
        L2_norms[quantity].append(
                # root_mean_square(
                # solution_family_array[quantity][i]
                # )
                max(solution_family_array[quantity][i])
                )
    for i, lambda_value in enumerate(lambda_value_array1):
        L2_norms1[quantity].append(
                # root_mean_square(
                # solution_family_array[quantity][i]
                # )
                max(solution_family_array1[quantity][i])
                )
    for i, lambda_value in enumerate(lambda_value_array2):
        L2_norms2[quantity].append(
                # root_mean_square(
                # solution_family_array[quantity][i]
                # )
                max(solution_family_array2[quantity][i])
                )

    axes_dict[quantity].plot(lambda_value_array, L2_norms[quantity], label="Hadamard point-splitting")
    axes_dict[quantity].plot(lambda_value_array1, L2_norms1[quantity], label="Mode sum formula")
    axes_dict[quantity].plot(lambda_value_array2, L2_norms2[quantity])

    axes_dict[quantity].legend(loc="best")

fontsize = 13
ax_rho.set_ylabel(r"max($\rho$)", fontsize=fontsize)
ax_rho.set_xlabel(r"$\lambda$", fontsize=fontsize)
ax_A0_induced.set_ylabel(r"max($-\int \int \rho$)", fontsize=fontsize)
# ax_rho.set_ylabel(r"$\sqrt{\langle \rho^2\rangle}$")
# ax_A0_induced.set_ylabel(r"$\sqrt{\langle A_0^2\rangle}$")



fig.tight_layout()
plt.show()
