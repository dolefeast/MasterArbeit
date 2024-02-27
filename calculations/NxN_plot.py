from scipy import special
from math_objects import system
from math_objects.unique_floats import float_in_array, unique_floats
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

LAMBDA_VALUE = 20
TOL = 1e-1
N_POINTS = 20
SCALAR_VALUE = 5
SCALAR_MASS = 9

n_images = 3
lower = 1
higher = 10
fig, ax = plt.subplots(n_images + 1, n_images + 1, figsize=(16, 9))

for lambda_index, lambda_value in enumerate(
    range(lower, higher, (higher - lower) // n_images)
):
    for scalar_index, scalar_mass in enumerate(
        range(lower, higher, (higher - lower) // n_images)
    ):
        print(f"Calculating for lambda={lambda_value}, mass={scalar_mass}")
        systesm = system.System(
            scalar_name="phi",
            n_points=N_POINTS,
            scalar_value=SCALAR_VALUE,
            field_strengh=lambda_value,
            scalar_mass=scalar_mass,
        )

        x = systesm.z
        omegas = []
        for omega_n in range(1, 100, 30):
            solution = sp.integrate.solve_bvp(
                systesm.differential_equation,
                systesm.dirichlet_boundary_conditions,
                systesm.z,
                (systesm.phi.value, systesm.phi.gradient.value),
                p=(omega_n,),
                verbose=0,
                max_nodes=3000,
                tol=TOL,
            )

            if not solution.success:
                fmt = "--"
            else:
                fmt = "-"

            if float_in_array(solution.p[0], omegas):
                continue
            else:
                ax[lambda_index, scalar_index].plot(
                    x,
                    (solution.p[0] - systesm.phi.charge * systesm.A0(x))
                    * np.real(solution.sol(x)[0]) ** 2
                    / scalar_mass
                    / lambda_value,
                    fmt,
                )
                ax[lambda_index, scalar_index].set_title(
                    f"lambda={lambda_value}, mass={scalar_mass}"
                )
                # plt.plot(x, solution.sol(x)[1])
                omegas.append(solution.p[0])


omegas = unique_floats(omegas, precision=1e-5)
print(sorted(omegas))
# ni = 1j*systesm.phi.mass**2/lambda_value/2 - 1/2
phi = systesm.phi

phi_n_dirichlet = lambda z, n: (phi.mass**2 + np.pi**2 * n**2) ** (-1 / 4) * (
    np.sin(np.pi * n * z)
    + lambda_value
    * np.sqrt(phi.mass**2 + np.pi**2 * n**2)
    / (2 * np.pi * np.abs(n))
    * (
        (1 / 2 - z) / np.pi / n * np.sin(np.pi * n * z)
        - z * (1 - z) * np.cos(np.pi * n * z)
    )
)

# for n in range(1, 3):
#    omega = np.sqrt(phi.mass**2 + np.pi**2 * n**2)
#    plt.plot(x, phi_n_dirichlet(x, n), label=f'Perturbative solution for n = {n}, omega = {omega}' )
#
# plt.legend()
plt.tight_layout()
plt.show()
