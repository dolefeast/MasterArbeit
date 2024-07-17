import matplotlib.pyplot as plt
from mpmath import quad
import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicSpline

def sign(x):
    # Returns sign of x
    return 2 * (x>0) - 1 

def float_in_array(x, array, tol=1e-3):
    # Checks if x is in given array up to some tolerance
    for element in array:
        if abs(x - element) < tol: 
            return True
    return False


def A0(z):
    # This is horrible I know. Usually all this is done inside a class
    # A0_induced should be a CubicSpline
    return - lambda_value * (z- 1/2) + A0_induced(z)

def eigenvalue_guess(n):
    # energy of a massless free particle
    return n * np.pi

def dirichlet_eigenstate(z, omega_n):
    # A first order approximation of a solution guess for the Klein-Gordon equation
    n = sign(omega_n) / np.pi

    if n == 0:
        return 0

    eigenstate = abs(omega_n) ** -(1 / 2) * (
        np.sin(np.pi * n * z)
        + lambda_value
        * omega_n
        / (2 * np.pi * abs(n))
        * ((1 / 2 - z) / np.pi / n * np.sin(np.pi * n * z) - z * (1 - z) * np.cos(np.pi * n * z))
    )

    return eigenstate

def dirichlet_eigenstate_gradient(z, omega_n):
    n = sign(omega_n) / np.pi

    if n == 0:
        return 0

    eigenstate_gradient = abs(omega_n) ** -(1 / 2) * (
        np.cos(np.pi * n * z) * np.pi * n
        + lambda_value
        * omega_n
        * (
            (
                -1 / np.pi / n * np.sin(np.pi * n * z)
                + (1 / 2 - z) * np.cos(np.pi * n * z)  # Deriv of first term
                - (1 - z) * np.cos(np.pi * n * z)
                + z * np.cos(np.pi * n * z)
            )
            + np.pi * n * z * (1 - z) * np.sin(np.pi * n * z)
        )
        / (2 * np.pi * abs(n))
    )

    return eigenstate_gradient

def Klein_Gordon_equation(
    z_array,
    f,
    p,
    ):
    # The differential equation
    omega_n = p[0] # Eigenvalue of the nth mode

    background_field = A0(z_array)

    klein_gordon = np.array(
        (f[1], -((omega_n - background_field) ** 2) * f[0])
    )

    return klein_gordon

def Dirichlet_bcs(ya, yb, p):
    # Since we are finding the eigenvalue,
    # need an additional BC to solve the equation
    bcs = np.array((ya[0], yb[0], yb[1] - 1))
    return bcs

def calculate_eigenstates(
        z,
        max_N,
        ):

    eigenvalues = []
    eigenstates = []
    eigenstate_gradients = []

    for n in range(-max_N, max_N+1):
        if n == 0:
            continue # Non physical case

        omega_n = eigenvalue_guess(n) # the eigenvalue guess

        # The guess for the solution, and its derivative, as
        # its solving a 2nd order differential eq
        y = (
                dirichlet_eigenstate(z, omega_n),
                dirichlet_eigenstate_gradient(z, omega_n),
                )
        eigenstate = solve_bvp(
                Klein_Gordon_equation,
                Dirichlet_bcs,
                z, # The x-axis,
                y,
                p=(omega_n,)
                )
        omega_n = eigenstate.p[0]

        if not eigenstate.success:
            # Sometimes the ODE doesnt converge
            raise Exception(f"Solution associated with omega_n={omega_n} did not converge")

        if float_in_array(omega_n, eigenvalues):
            # Sometimes an eigenvalue converges to another possible solution
            raise Exception(f"Found repeated eigenvalue {omega_n}")

        eigenvalues.append(omega_n)
        eigenstates.append(eigenstate.sol(z)[0])
        eigenstate_gradients.append(eigenstate.sol(z)[1])

    return {
            "eigenvalues": eigenvalues,
            "eigenstates": eigenstates,
            "eigenstate_gradients": eigenstate_gradients,
            }

def normalize_eigenstate_family(z, solution_family:dict):
    # Due to the boundary conditions, the solutions are not normalized
    # This assumes that the A0 is somewhere defined.

    # Warning, math
    for n, (eigenvalue, eigenstate, eigenstate_gradient) in enumerate(
            zip(
            solution_family["eigenvalues"],
            solution_family["eigenstates"],
            solution_family["eigenstate_gradients"],
            )
            ):

        eigenstate = CubicSpline(z, eigenstate)

        def rho_n_without_normalizing(z):
            # Normalizing wrt the symplectic norm
            # the solutions need not be real.
            return (eigenvalue - A0(z)) * abs(eigenstate(z))**2

        # Calculate the norm
        norm_squared = abs(float(quad(rho_n_without_normalizing, [0, 1])))

        solution_family["eigenstates"][n] = eigenstate(z)/np.sqrt(norm_squared)
        solution_family["eigenstate_gradients"][n] = eigenstate_gradient/np.sqrt(norm_squared)

    return solution_family

def calculate_rho(solution_family:dict):
    eigenvalues = solution_family["eigenvalues"]
    eigenstates = solution_family["eigenstates"]


def convolution(rho, max_N, n_points):
    # This depends on pretty hardcore fine tuning.
    # As long as n_points is exactly 8 * (max_N + 1),
    # the filtering should work fine
    delta_N = n_points // ( max_N + 1) + 1
    convoluting_array = [0] * delta_N

    # fine tuned parameters
    surround_peaks = 0.00
    peaks = 0.4
    middle = 6
    edges = 3

    # On the borders
    array[0] = edges
    array[-1] = edges

    # Peak of the sine
    array[delta_N // 4 + 1] = surround_peaks
    array[delta_N // 4 - 1] = surround_peaks

    # Trough of the sine
    array[3 * (delta_N // 4) + 1] = surround_peaks
    array[3 * (delta_N // 4) - 1] = surround_peaks

    array[delta_N // 4] = peaks
    array[3 * (delta_N // 4)] = peaks

    # In the middle of the curve
    array[delta_N // 2] = middle

    # Normalizing it so that we don't get extra factors
    array = np.array(array) / sum(array)

    rho_filtered = np.convolve(rho, array, mode="same")
 
    return rho_filtered # This will stield yield some noise as 
                      # this is not the full smoothing algorithm.


fig, ax = plt.subplots()

n_points = 200
z = np.linspace(0, 1, n_points)
max_N = 24
lambda_value = 2
A0_induced = CubicSpline(z, np.zeros_like(z))

solution_family = calculate_eigenstates(z, max_N)

solution_family = normalize_eigenstate_family(z, solution_family)

for n in range(2*max_N-1):
    ax.plot(solution_family["eigenstates"][n], label=f"n = {n}")
    ax.legend()
    fig.canvas.draw()
    plt.pause(0.001)
    plt.cla()
    input("aa")


# solution_family = normalize_eigenstate_family(solution_family)
# plt.plot(solution_family["eigenstates"][-1](z))

# plt.show()
