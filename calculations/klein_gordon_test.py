from scipy import special
from math_objects import system
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

lambda_value = 20
TOL = 1e-2
N_POINTS = 100
scalar_value = 1
scalar_mass = 20

fig, ax = plt.subplots(figsize=(16, 9))

systesm = system.System(scalar_name = 'phi', 
        n_points = N_POINTS,
        scalar_value = scalar_value,
        field_strengh = lambda_value,
        scalar_mass = scalar_mass)

eigenstates = []
omega_solution_array = []
x = systesm.z

total_charge_density = 0
for omega_n in range(0,150,6):
    omega_n /= 1
    solution = sp.integrate.solve_bvp(systesm.differential_equation,
            systesm.dirichlet_boundary_conditions,
            systesm.z, 
            (systesm.phi.value, systesm.phi.gradient.value),
            p=(omega_n, ),
            verbose=0,
            max_nodes=3000,
            tol=TOL)

    if not solution.success:
        fmt = '--'
    else:
        fmt = '-'

    #solution_field = normalize(lambda z: solution.sol(z)[0])
    solution_squared = lambda z: np.real(solution.sol(z)[0]**2)
    solution_norm_squared = sp.integrate.quad(solution_squared, 0, 1)[0]

    charge_density_n = (solution.p[0]- systesm.phi.charge*systesm.A0(x))*np.real(solution.sol(x)[0])**2/scalar_mass/lambda_value/solution_norm_squared
    total_charge_density += charge_density_n

    ax.plot(x,  charge_density_n, fmt)
    ax.set_title(f'lambda={lambda_value}, mass={scalar_mass}')
    omega_solution = solution.p[0]
    eigenstates.append(solution)
    omega_solution_array.append(omega_solution)


#plt.plot(x, systesm.calculate_charge_density(x, eigenstates))


ax.plot(x,  total_charge_density, '-')
plt.tight_layout()
plt.show()
