from math_objects.system import System
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
from math_objects.savitzky_golay import savitzky_golay

import scipy as sp
from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt

lambda_value = 1
TOL = 1e-2
N_POINTS = 1000
scalar_value = 1
scalar_mass = 1
scalar_charge = 1

fig, (ax_fields, ax_omegas) = plt.subplots(1, 2, figsize=(16, 9))

system = System(scalar_name = 'phi', 
        n_points = N_POINTS,
        scalar_charge = scalar_charge,
        scalar_value = scalar_value,
        field_strength = lambda_value,
        scalar_mass = scalar_mass)

#system.phi.value = np.sin(system.z)

eigenstates = []
omega_solution_array = []

total_charge_density = 0

omega_boundary = 550
n_of_solutions = int(omega_boundary/3)
eigenstate_array = system.calculate_N_eigenstates(n_of_solutions, 
       -omega_boundary, 
        omega_boundary)


total_charge_density = system.calculate_total_charge_density(eigenstate_array)

for i, eigenstate in enumerate(eigenstate_array):
    if not eigenstate.success:
        fmt = 'r'
    else:
        fmt = 'g'

    omega_solution = eigenstate.p[0]
    if float_in_array(omega_solution, omega_solution_array):
        continue
    else:

        #ax_fields.plot(system.z,  system.calculate_charge_density(eigenstate)(system.z), fmt, alpha=0.5, linewidth=0.6)
        ax_omegas.scatter(i,  omega_solution, color=fmt)
        ax_omegas.set_title(f'$\lambda={lambda_value}, m={scalar_mass}$')
        omega_solution_array.append(omega_solution)

total_charge_density_array = total_charge_density(system.z)

ax_fields.plot(
        system.z,
        total_charge_density_array,
        label='Rough charge density'
        )

#SMOOTHING
signal = savitzky_golay(
        total_charge_density_array,
        5*system.n_points/(n_of_solutions+1), 
        1
        )
ax_fields.plot(
        system.z,
        signal,
        label='Savitzky golay smoothing'
        )

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

signal = smooth(
        total_charge_density_array,
        n_of_solutions//5
        )

ax_fields.plot(system.z, signal, label='convolution smoothing')
ax_fields.legend(loc='best')

#electric_field = system.new_electric_field(eigenstate_array)
#vector_field = system.new_vector_field(eigenstate_array)
#print(vector_field)
#ax_fields.plot(vector_field.t, vector_field.y[0], label='EM field')


#ax_fields.plot(system.z,  total_charge_density, '-', linewidth=2)

plt.tight_layout()
plt.show()
