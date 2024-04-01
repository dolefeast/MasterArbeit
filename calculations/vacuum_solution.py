from math_objects import Vacuum_Solution
from math_objects.savitzky_golay import savitzky_golay

import numpy as np
import matplotlib.pyplot as plt

lambda_value = 0.6
TOL = 1e-2
N_POINTS = 1000
m = 1
e = 1

fig, (ax_fields, ax_omegas) = plt.subplots(1, 2, figsize=(16, 9))

system = Vacuum_Solution(
    lambda_value=lambda_value,
    scalar_name="phi",
    n_points=N_POINTS,
    e=e,
    m=m,
)

omega_boundary = 150
n_of_solutions = omega_boundary * 3

# system.phi.value = np.sin(system.z)

eigenstates = []
omega_solution_array = []

eigenstate_array = system.calculate_eigenstates()

total_charge_density_callable = system.calculate_total_charge_density(eigenstate_array)

for i, eigenstate in enumerate(eigenstate_array):
    if not eigenstate.success:
        fmt = "r"
    else:
        fmt = "g"

    omega_solution = eigenstate.p[0]

    # ax_fields.plot(system.z,  system.calculate_charge_density(eigenstate)(system.z), fmt, alpha=0.5, linewidth=0.6)
    ax_fields.plot(
        system.z,
        eigenstate.sol(system.z)[0] / omega_solution,
        fmt,
        alpha=0.5,
        linewidth=0.6,
    )
    ax_omegas.scatter(
        np.sqrt(omega_solution**2 - m**2) / np.pi, omega_solution, marker="x"
    )
    ax_omegas.set_title(f"$\lambda={lambda_value}, m={m}$")
    omega_solution_array.append(omega_solution)

total_charge_density_array = total_charge_density_callable(system.z)
total_charge_density_array[0] = 0

ax_fields.plot(
    system.z, total_charge_density_array, label="Rough charge density", alpha=0.6
)

# SMOOTHING
smoothed_charge_density = savitzky_golay(
    total_charge_density_array, 10 * system.n_points // (n_of_solutions + 1), 2
)

ax_fields.plot(
    system.z, smoothed_charge_density, label="Savitzky golay smoothing 1st time"
)

smoothed_charge_density = savitzky_golay(
    smoothed_charge_density, 10 * system.n_points // (n_of_solutions + 1), 4
)
ax_fields.plot(
    system.z, smoothed_charge_density, label="Savitzky golay smoothing 2nd time"
)
ax_fields.plot(system.z, system.A0(system.z), label="electric potential")

# electric_field = system.new_electric_field(eigenstate_array)
# vector_field = system.new_vector_field(smoothed_charge_density)
##ax_fields.plot(system.z, vector_field.y[0], label='$A_0$ field')
# ax_fields.plot(system.z, vector_field.y[0], label='$A_0$ field')

ax_fields.legend(loc="best")

plt.tight_layout()
plt.show()
