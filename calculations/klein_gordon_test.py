from math_objects.system import System
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
from math_objects.savitzky_golay import savitzky_golay

import scipy as sp
from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt

lambda_value = 10
TOL = 1e-2
N_POINTS = 1000
scalar_value = 1
scalar_mass = 5
scalar_charge = 1

fig, (ax_fields, ax_omegas) = plt.subplots(1, 2, figsize=(16, 9))

system = System(
    scalar_name="phi",
    n_points=N_POINTS,
    scalar_charge=scalar_charge,
    scalar_value=scalar_value,
    field_strength=lambda_value,
    scalar_mass=scalar_mass,
)

# system.phi.value = np.sin(system.z)

eigenstates = []
omega_solution_array = []

omega_boundary = 500
n_of_solutions = omega_boundary // 2
eigenstate_array = system.calculate_N_eigenstates(
    n_of_solutions, -omega_boundary, omega_boundary
)


total_charge_density_callable = system.calculate_total_charge_density(eigenstate_array)

for i, eigenstate in enumerate(eigenstate_array):
    if not eigenstate.success:
        fmt = "r"
    else:
        fmt = "g"

    omega_solution = eigenstate.p[0]
    if float_in_array(omega_solution, omega_solution_array):
        continue

    # ax_fields.plot(system.z,  system.calculate_charge_density(eigenstate)(system.z), fmt, alpha=0.5, linewidth=0.6)
    ax_omegas.scatter(i, omega_solution, color=fmt)
    ax_omegas.set_title(f"$\lambda={lambda_value}, m={scalar_mass}$")
    omega_solution_array.append(omega_solution)

total_charge_density_array = total_charge_density_callable(system.z)
total_charge_density_array[0] = 0

ax_fields.plot(system.z, total_charge_density_array, label="Rough charge density", alpha=0.6)

# SMOOTHING
print("Calculating forwards smoothing")
signal = savitzky_golay(
    total_charge_density_array, 10*system.n_points // (n_of_solutions + 1), 2
)

print("Calculating backwards smoothing")
signal_backwards = savitzky_golay(
        total_charge_density_array[::-1], 10*system.n_points // (n_of_solutions + 1), 0
)

ax_fields.plot(system.z, signal, label="Savitzky golay smoothing 1st time")
ax_fields.plot(system.z, signal_backwards, label="Savitzky golay backwards smoothing 1st time")

signal = savitzky_golay(
    signal, 10*system.n_points // (n_of_solutions + 1), 4
)
# ax_fields.plot(system.z, signal, label="Savitzky golay smoothing 2nd time")

#box_convolution(total_charge_density_array)
ax_fields.legend(loc="best")

# electric_field = system.new_electric_field(eigenstate_array)
#vector_field = system.new_vector_field(total_charge_density_callable)
#print(vector_field)
#ax_fields.plot(vector_field.t, vector_field.y[0], label='EM field')

plt.tight_layout()
plt.show()
