from math_objects import Vacuum_Solution
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
from math_objects.savitzky_golay import savitzky_golay

import scipy as sp
from scipy.signal import savgol_filter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'full') / w

N_mode_cutoff = 150

lambda_value = 1
TOL = 1e-2
N_POINTS = 1000
m = 4.7
e = 1

fig, (ax_fields, ax_omegas) = plt.subplots(1, 2, figsize=(16, 9))
fig.suptitle(f"$\lambda={lambda_value}, m={m}$")

system = Vacuum_Solution(
    lambda_value=lambda_value,
    scalar_name="phi",
    n_points=N_POINTS,
    N_mode_cutoff=N_mode_cutoff+1,
    e=e,
    m=m,
    float_tol=1e-1
)

# system.phi.value = np.sin(system.z)

eigenstates = []
omega_solution_array = []

eigenstate_array = system.calculate_eigenstates()

z = system.z # Because I consistently make the mistake of only asking for z

for i, eigenstate in enumerate(eigenstate_array):
    omega_solution = eigenstate.p[0]
    if omega_solution < 0:
        omega_scatter_color = 'r'
        omega_scatter_marker = 'x'
    elif omega_solution > 0:
        omega_scatter_color = 'g'
        omega_scatter_marker = 'o'

    ax_omegas.scatter(
             np.sqrt(omega_solution**2 - m**2)/np.pi + 0.3*(omega_solution<0),
             np.abs(omega_solution),
             marker=omega_scatter_marker,
             color=omega_scatter_color
             )
 
    omega_solution_array.append(omega_solution)

total_charge_density = system.calculate_total_charge_density(
        eigenstate_array, 
            renormalization=False,
            smoothing=False,
            )(z) # calculate_... returns a function

ax_fields.plot(
        total_charge_density,
        label='Unfiltered charge density'
        )

window_initial = 50
window_max = 600
window_frequency = 5
ax_fields.plot(
       moving_average(total_charge_density, window_initial),
       label=f'Moving average filtering. window={window_initial}'
       )

while True:
    plt.draw()
    while True:
        try:
            val = int(input('Input desired window size: '))
            break
        except ValueError:
            print('Input must be an integer')
    ax_fields.clear()
    ax_fields.plot(
            moving_average(
                moving_average(
                    total_charge_density,
                    int(val)//5+1
                    ),
                int(val),
                ),
            label=f'Moving average filtering. window={window_initial}'
            )
    plt.pause(0.01)

ax_fields.set_xlabel('z')
ax_fields.set_ylabel(r'$\rho(z)$')

ax_omegas.set_xlabel('n')
ax_omegas.set_ylabel(r'$\omega_n$')

ax_omegas.plot([],[], 'go', label='positivw eigenvalue')
ax_omegas.plot([],[], 'rx', label='negative eigenvalue')
ax_fields.legend(loc='best')
ax_omegas.legend(loc='best')
plt.tight_layout()
plt.show()
