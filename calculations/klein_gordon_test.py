from scipy import special
from math_objects import system
from math_objects.unique_floats import float_in_array, unique_floats
from math_objects.normalize import normalize
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

lambda_value = 5
TOL = 1e-2
N_POINTS = 300
scalar_value = 1
scalar_mass = 1

fig, (ax_fields, ax_omegas) = plt.subplots(1, 2, figsize=(16, 9))

systesm = system.System(scalar_name = 'phi', 
        n_points = N_POINTS,
        scalar_value = scalar_value,
        field_strengh = lambda_value,
        scalar_mass = scalar_mass)

eigenstates = []
omega_solution_array = []
x = systesm.z

total_charge_density = 0

omega_limit = 50

for i, omega_guess in enumerate(np.linspace(-omega_limit, omega_limit, 500)):
    solution = sp.integrate.solve_bvp(systesm.differential_equation,
            systesm.dirichlet_boundary_conditions,
            systesm.z, 
            (systesm.phi.value, systesm.phi.gradient.value),
            p=(omega_guess, ),
            verbose=0,
            max_nodes=3000,
            tol=TOL)

    if not solution.success:
        fmt = 'r'
    else:
        fmt = 'g'

    #solution_field = normalize(lambda z: solution.sol(z)[0])
    solution_squared = lambda z: (solution.p[0] - systesm.phi.charge * systesm.A0(z))*np.real(solution.sol(z)[0]**2)
    solution_norm_squared = sp.integrate.quad(solution_squared, 0, 1)[0]

    #print(f'Guessing for omega = {omega_n}\n\tomega_n = {solution.p[0]} \n\tsolution_norm_squared={solution_norm_squared}')
    ax_omegas.scatter(i,  omega_guess, color='b')
    omega_solution = solution.p[0]
    if float_in_array(omega_solution, omega_solution_array):
        continue
    else:
        field_factor = np.sign(omega_solution)*(omega_solution- systesm.phi.charge*systesm.A0(x))
        charge_density_n = field_factor*np.real(solution.sol(x)[0])**2/scalar_mass/lambda_value/solution_norm_squared
        total_charge_density += charge_density_n

        ax_fields.plot(x,  charge_density_n, fmt, alpha=0.5, linewidth=0.6)
        ax_omegas.scatter(i,  omega_solution, color=fmt)
        ax_omegas.set_title(f'$\lambda={lambda_value}, m={scalar_mass}$')
        eigenstates.append(solution) #Watch out there are repeating eigenstates.
        omega_solution_array.append(omega_solution)

    


#plt.plot(x, systesm.calculate_charge_density(x, eigenstates))


ax_fields.plot(x,  total_charge_density, '-', linewidth=2)
plt.tight_layout()
plt.show()
