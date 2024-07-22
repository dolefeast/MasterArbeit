import matplotlib.pyplot as plt
from mpmath import quad
import numpy as np
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicSpline

from itertools import cycle

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
    return - lambda_value * (z- 1/2) + a * A0_induced(z)

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
    # need an additional BC to solve the equation.
    # y'(1) = 1 is arbitrary but preserves generality
    # as long as the solution to the ODE is later normalized
    bcs = np.array((ya[0], yb[0], yb[1] - 1))
    return bcs

def calculate_eigenstates(max_N):

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
                p=(omega_n,),
                max_nodes=9e19,
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

def normalize_eigenstate_family(solution_family:dict):
    # Due to the boundary conditions, the solutions are not normalized
    # This assumes that the A0 is somewhere defined.

    # Had to initiate the dict again, I am not sure why
    # If I don't do it, it doesn't really update.
    eigenvalues = solution_family["eigenvalues"]
    eigenstates = []
    eigenstate_gradients = []
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

        norm = np.sqrt(norm_squared)

        eigenstates.append(eigenstate(z)/norm)
        eigenstate_gradients.append(eigenstate_gradient/norm)

    return {
            "eigenvalues":eigenvalues,
            "eigenstates":eigenstates,
            "eigenstate_gradients":eigenstate_gradients,
            }

def calculate_rho(solution_family:dict):
    # If the solutions are not normalized this does not work
    # Again, this assumes that A0 is defined somewhere before this function is called
    eigenvalues = solution_family["eigenvalues"]
    eigenstates = solution_family["eigenstates"]

    # The total charge density. 
    # Initiate as np.array to add them
    rho = np.zeros_like(z) 

    for n, (eigenvalue, eigenstate) in enumerate(
            zip(
                eigenvalues,
                eigenstates,
                )
            ):
        # The charge density associated to the nth mode
        # eigenstate is np array. A0 is a callable. 
        # A0(z) is of the same shape as eigenstate
        rho_n = (eigenvalue - A0(z)) * abs(eigenstate) ** 2
        rho += rho_n

    return rho


def filter_rho(rho):
    """
    Parameters:
        rho: a noisy signal (with a very particular kind of noise)
    Returns:
        rho_filtered: the filtered signal.
    """
    # This exploits the fact that we know exactly what the 
    # noise looks like. It can be convoluted away with
    # a very specific function. Not particularly
    # relevant.

    # As long as n_points is exactly 8 * (max_N + 1),
    # the filtering should work fine
    # As with A0, n_points, max_N are necessarily
    # defined somewhere before the function is called.
    delta_N = n_points // ( max_N + 1) + 1
    convoluting_array = [0] * delta_N

    # Some fine tuned parameters
    surround_peaks = 0.00
    peaks = 0.4
    middle = 6
    edges = 3

    # On the borders
    convoluting_array[0] = edges
    convoluting_array[-1] = edges

    # Peak of the sine
    convoluting_array[delta_N // 4 + 1] = surround_peaks
    convoluting_array[delta_N // 4 - 1] = surround_peaks

    # Trough of the sine
    convoluting_array[3 * (delta_N // 4) + 1] = surround_peaks
    convoluting_array[3 * (delta_N // 4) - 1] = surround_peaks

    convoluting_array[delta_N // 4] = peaks
    convoluting_array[3 * (delta_N // 4)] = peaks

    # In the middle of the curve
    convoluting_array[delta_N // 2] = middle

    # Normalizing it so that we don't get extra factors
    convoluting_array = np.array(convoluting_array) / sum(convoluting_array)

    rho_filtered = np.convolve(rho, convoluting_array, mode="same")
 
    return rho_filtered # This will stield yield some noise as 
                      # this is not the full smoothing algorithm.

def calculate_A0_induced(rho):
    rho_interpolated = CubicSpline(z, rho)
    # Math: 
    # Since gauge freedom, I can choose the offset of the induced field.
    # Due to the antisymmetry of the induced rho, I can choose this
    # offset so that the A0 field is antisymmetric with respect to
    # the midpoint z=1/2
    A0_induced_shifted = rho_interpolated.antiderivative(2)
    offset = A0_induced_shifted(1/2)
    # This forces A0_induced(1/2) = 0
    # The minus sign is because physics.
    A0_induced = lambda z: -(A0_induced_shifted(z) - offset)

    return A0_induced

def plot_rho_A0_induced(
        rho, 
        A0_induced,
        fmt,
        color,
        alpha,
        ):
    ax_rho.plot(z, rho, fmt, color=color, alpha=alpha)
    ax_A0_induced.plot(z, A0_induced, fmt, color=color, alpha=alpha)


def calculate_A0_total(A0_induced):
    # A0_induced is a function
    return lambda z: -lambda_value * (z - 1/2) + a*A0_induced(z)


if __name__ == "__main__":

    plot = True
    if plot:
        fig, (ax_rho, ax_A0_induced, ax_eigenvalues) = plt.subplots(3)
        eigenvalues_array = []
        lambda_array = []

    # Number of modes
    n_points = 200
    # x-axis
    z = np.linspace(0, 1, n_points)

    a = 1
    max_N = 1 #24 # Cutoff for the mode solutions
               # Note that there will be 2 * max_N - 1 solutions
    lambda_min = 7 # Strength of the external electric field
    lambda_max = 25 # Strength of the external electric field
    lambda_step=0.25

    n_iterations = 15 # Number of iterations for each value of lambda

    # Need to initialize the A0_induced
    A0_induced = CubicSpline(z, np.zeros_like(z))
    color_cycle = cycle(["#FE5F55","#3423A6","#710627","#7180B9","#163B20"])

    for i, lambda_value in enumerate(np.arange(lambda_min, lambda_max, lambda_step)):
        print('λ =', lambda_value)
        color = next(color_cycle) # To distinguish between different lambda values
        if lambda_value == 17.0:
            print('We got here, increasing n_iterations')
            n_iterations *= 2
        try:
            for n in range(n_iterations):
                solution_family = calculate_eigenstates(max_N)

                normalized_solution_family = normalize_eigenstate_family(solution_family)

                rho = calculate_rho(normalized_solution_family)
                # rho = filter_rho(rho)
                A0_induced = calculate_A0_induced(rho)

                A0 = calculate_A0_total(A0_induced)
                if plot and not i%3:
                    plot_rho_A0_induced(rho, A0_induced(z), '--', color=color, alpha=0.5)
                if not n: # For the first iteration
                    print('n=0, max(A0_induced) =', max(A0_induced(z)))
                
            if plot:
                eigenvalues_array.append(solution_family["eigenvalues"])
                lambda_array.append(lambda_value)
                plot_rho_A0_induced(rho, A0_induced(z), '', color=color, alpha=1)
            # print(f'n={n_iterations}, max(A0_induced) =', max(A0_induced(z)))
            A0 = calculate_A0_total(lambda z: (lambda_value + lambda_step)/lambda_value * A0_induced(z))

        except Exception as exception:
            e = exception
            break
        except KeyboardInterrupt as exception:
            e = exception
            break

    if plot:
        ax_eigenvalues.plot(lambda_array, eigenvalues_array, 'b')
        plt.show()
    
    raise e

