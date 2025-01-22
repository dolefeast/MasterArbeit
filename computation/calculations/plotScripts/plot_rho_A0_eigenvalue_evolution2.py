from read_data_scripts import get_Posix_for_quantities, open_Posix_dict
from read_files import get_directory_m_a, str_to_float

from itertools import cycle
import sys
import matplotlib.pyplot as plt
import numpy as np

def get_index_for_value_from_below(array, value):
    """
    Given an array, return the biggest index i for which array[i] <= value
    """
    for i, element in enumerate(array):
        if element > value:
            return i-1

try:
    filter_regex = sys.argv[1]
except IndexError:
    filter_regex = ""


bcs = "dirichlet"
max_lambda_density = None

directory, m, a = get_directory_m_a(bcs, filter_regex=filter_regex)
directory1, m1, a1 = get_directory_m_a(bcs, filter_regex=filter_regex)
directory2, m2, a2 = get_directory_m_a(bcs, filter_regex=filter_regex)

desired_quantities = ["rho", "A0_induced"]
posix_dict = get_Posix_for_quantities(m, a, directory=directory, max_lambda_density=max_lambda_density)
solution_family_array = open_Posix_dict(posix_dict, desired_quantities=["eigenvalues"] + desired_quantities)

posix_dict1 = get_Posix_for_quantities(m, a, directory=directory1, max_lambda_density=max_lambda_density)
solution_family_array1 = open_Posix_dict(posix_dict1, desired_quantities=["eigenvalues"] + desired_quantities)

posix_dict2 = get_Posix_for_quantities(m2, a2, directory=directory2, max_lambda_density=max_lambda_density)
solution_family_array2 = open_Posix_dict(posix_dict2, desired_quantities=["eigenvalues"] + desired_quantities)


maxN = len(solution_family_array["eigenvalues"][0])//2
if maxN == 1: 
    omega_ylim = 1.2*solution_family_array["eigenvalues"][0][:]
    omega_ylim[0] = 0
else:
    omega_ylim = (0, 10)

lambda_value_array = solution_family_array["lambda_value"]
lambda_value_array1 = solution_family_array1["lambda_value"]
lambda_value_array2 = solution_family_array2["lambda_value"]

print("The maximum lambda value obtained is", lambda_value_array[-1]) #, ", with a minimum lambda step of", lambda_value_array[-1]-lambda_value_array[-2])

fig1 = plt.figure(f'{directory}, m={str_to_float(1, m)}, a={str_to_float(1, a)}1')
fig2 = plt.figure(f'{directory}, m={str_to_float(1, m)}, a={str_to_float(1, a)}2')
fig3 = plt.figure(f'{directory}, m={str_to_float(1, m)}, a={str_to_float(1, a)}3')

ax_rho = fig1.subplots()
ax_A0_induced = fig2.subplots()
ax_omega = fig3.subplots()

axes_dict = {"rho":ax_rho, "A0_induced":ax_A0_induced}

n_solutions = len(lambda_value_array )

interesting_lambda_value_array = [5]

colors = cycle(["#0075F2",
        "#CD4631",
        "#9A348E",])


for quantity in desired_quantities:
    axes_dict[quantity].set_xlabel("z")
    for interesting_lambda_value in interesting_lambda_value_array:
        index = get_index_for_value_from_below(lambda_value_array, interesting_lambda_value)
        index1 = get_index_for_value_from_below(lambda_value_array1, interesting_lambda_value)

        recovered_quantity = solution_family_array[quantity][index]
        recovered_quantity1 = solution_family_array1[quantity][index1]

        n_points = len(recovered_quantity)
        n_points1 = len(recovered_quantity1)
        z = np.linspace(0, 1, n_points) 
        z1 = np.linspace(0, 1, n_points1) 

        # alpha = 0.2 + 0.8 * ((i+1)/len(lambda_value_array))**3
        color = next(colors)

        axes_dict[quantity].plot(
                z,
                recovered_quantity,
                # 'b',
                # alpha=alpha,
                label=r"Mode sum formula", 
                )
        axes_dict[quantity].plot(
                z1,
                recovered_quantity1,
                # 'b',
                # alpha=alpha,
                label=f"Hadamard point-split renormalization", 
                )

ax_omega.plot(lambda_value_array, [omega_array[1] for omega_array in solution_family_array["eigenvalues"]], label=r"$\omega_1^\lambda$ as calculated in [AW83]")
ax_omega.plot(lambda_value_array1, [omega_array[12] for omega_array in solution_family_array1["eigenvalues"]], label=r"$\omega_1^\lambda$ as calculated in [WZ20]")
ax_omega.plot(lambda_value_array2, [omega_array[1] for omega_array in solution_family_array2["eigenvalues"]], label=r"$\omega_1^\lambda$ without backreaction")

ax_rho.set_ylabel(r"$\rho$")
ax_A0_induced.set_ylabel(r"$A_0$")

ax_rho.set_title(r"$\lambda=$" + str(5))
ax_A0_induced.set_title(r"$\lambda=$" + str(5))

ax_rho.legend(loc="best")
ax_A0_induced.legend(loc="best")
ax_omega.legend(loc="best")

ax_omega.set_xlabel(r"$\lambda$")
ax_omega.set_ylabel(r"$\omega_n$")

ax_omega.set_ylim(omega_ylim)

fig1.tight_layout()
fig2.tight_layout()
fig3.tight_layout()

plt.show()
