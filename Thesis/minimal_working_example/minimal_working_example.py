import matplotlib.pyplot as plt
from mpmath import quad
from itertools import count
import numpy as np
from numpy import savetxt
from scipy.integrate import solve_bvp
from scipy.interpolate import CubicSpline
import os

from itertools import cycle

#TODO: 
# - Organize script into a single function

def sign(x):
    # Returns sign of x
    return 2 * (x>0) - 1 

def float_in_array(x, array, tol=1e-3):
    # Checks if x is in given array up to some tolerance
    for element in array:
        if abs(x - element) < tol: 
            return True
    return False

def float_to_str(value, sig_digs=3):
    if isinstance(value, str):
        return value
    return str(
        round(
            float(
                value
            ),
            sig_digs
        )
    ).replace(".", "_")

def str_to_float(value:str):
    return float(
        value.replace("_", ".")
        )

def root_mean_square(y):
    return np.sqrt(np.mean(y**2))

def A0(z):
    # This is horrible I know. Usually all this is done inside a class
    # A0_induced should be a CubicSpline
    return - lambda_value * (z- 1/2) + a * A0_induced(z)

def eigenvalue_guess(n):
    # energy of a massless free particle
    return n * np.pi

def dirichlet_eigenstate(z, omega_n):
    # A first order approximation of a solution guess for the Klein-Gordon equation
    if omega_n==0:
        return np.zeros_like(z)
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
    if omega_n==0:
        return np.zeros_like(z)
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

def fun_jac(z_array, y, p):
    omega = p[0]
    u, v = y # u = phi, v = phi'
    
    df_dy = np.zeros((len(y), len(y), len(z_array)))

    # Calculated by hand
    df_dy[0, 0, :] = np.zeros_like(z_array)
    df_dy[0, 1, :] = -(omega - A0(z_array))**2 + m**2
    df_dy[1, 0, :] = np.ones_like(z_array)
    df_dy[1, 1, :] = np.zeros_like(z_array)

    df_dp = np.zeros((len(y), len(p), len(z_array)))
    df_dp[0, 0, :] = np.zeros_like(z_array)
    df_dp[1, 0, :] = - 2 * omega * ( omega - e * A0(z_array) ) * u

    return df_dy, df_dp


def calculate_eigenstates(max_N,
        max_nodes=1e20,
        tol=1e-2,
        ):

    eigenvalues = []
    eigenstates = []
    eigenstate_gradients = []

    eigenvalue_guess_array = solution_family["eigenvalues"]
    eigenstate_guess_array = solution_family["eigenstates"]
    eigenstate_gradient_guess_array = solution_family["eigenstate_gradients"]

    for (
            eigenvalue_guess,
            eigenstate_guess,
            eigenstate_gradient_guess
            ) in zip(
                    eigenvalue_guess_array,
                    eigenstate_guess_array,
                    eigenstate_gradient_guess_array,
                    ):
        if eigenvalue_guess == 0:
            continue # Non physical case

        # omega_n = eigenvalue_guess[n] # the eigenvalue guess

        # The guess for the solution, and its derivative, as
        # its solving a 2nd order differential eq
        y = (
                eigenstate_guess,
                eigenstate_gradient_guess,
                )
        eigenstate = solve_bvp(
                Klein_Gordon_equation,
                Dirichlet_bcs,
                z, # The x-axis,
                y,
                p=(eigenvalue_guess,),
                max_nodes=max_nodes,
                fun_jac=fun_jac,
                tol=tol,
                )
        omega_n = eigenstate.p[0]

        if not eigenstate.success:
            # Sometimes the ODE doesnt converge
            #raise Exception(f"Solution associated with omega_n={omega_n} did not converge")
            raise Exception(eigenstate.message)

        if float_in_array(omega_n, eigenvalues) or abs(omega_n - eigenvalue_guess ) > 2*abs(eigenvalue_guess):
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
        rho += 1/2*rho_n

    return rho


def extend_signal(
        x,
        y,
        padding_size=None
        ):
    assert len(x) == len(y)

    # if padding_size is not specified,
    # then the padding is y itself
    if padding_size is None:
        padding_size = len(y)
    # wrap gives periodization of the function
    y_padded = np.pad(y, padding_size, mode="wrap")

    # need to extend x, too
    dx = x[1] - x[0]

    x_padded = np.linspace(-dx * padding_size, 1 + dx * padding_size, len(y_padded))

    return x_padded, y_padded

def remove_neighbourhood(
    x,
    y,
    points: (float)=(0,1),
    neighbourhood_size: float=None,
    force_zero: bool=True,
):
    """
    Given curve (x, y) with problematic neighbourhoods around points=(x1, x2, ...), take their neighbourhood with neighbourhood_size=neighbourhood_size away and interpolate around it, thus smoothing the curve.
    """

    if neighbourhood_size is None:
        # Take out the points corresponding to the convolution array
        neighbourhood_size = (
        1 / ( max_N + 1 ) 
        + 1 / n_points # Without this, the algorithm takes away
                            # an "open" neighbourhood (can't be open since it's discrete), 
                            # and adding it converts it into a closed one.
                            # Corresponds to adding a dz.
            )

    x_list = list(x)
    y_list = list(y)
    
    x_array = np.array(x)
    
    # A python array to use its addition properties
    idx = []

    # The points to take out of the array
    for p in points:
        
        window = np.where(abs(x_array - p) <= neighbourhood_size)[0].tolist()
        idx += window

    idx = np.reshape(
        np.array(idx),
        -1,  # 1-d array
    )
    x_list = [x for i, x in enumerate(x_list) if i not in idx]
    y_list = [y for i, y in enumerate(y_list) if i not in idx]

    if force_zero:  # Force the removed values to go through 0
        for p in points:
            for x_index, x_value in enumerate(x_list):
                if x_value > p:
                    x_list = x_list[:x_index] + [p] + x_list[x_index:]
                    y_list = y_list[:x_index] + [0] + y_list[x_index:]
                    break

    return np.array(x_list), np.array(y_list)

def remove_and_interpolate(
    x: [float],
    y: [float],
    points=(0, 1),
    neighbourhood_size:float=None,
    force_zero: bool=True,
):
    x_removed, y_removed = remove_neighbourhood(
        x,
        y,
        points=points,
        neighbourhood_size=neighbourhood_size,
        force_zero=True
    )

    interpolated_curve = CubicSpline(
            x_removed,
            y_removed,
            )

    return x, interpolated_curve(x)  # So that both are arrays

def return_to_0_1(
        x, 
        y):

    idx = np.where(
        np.logical_and(
            x >= 0,
            x <= 1,
        )
    )
    return x[idx], y[idx]

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
def extend_and_filter(
    x,
    y,
    neighbourhood_size=None,
    padding_size=None,
    points=(0, 1),
    force_zero=True,
):
    """
    Parameters:
    x: [float], the x-array of the signal. increasing and goes from 0 to 1
    y: [float], the y-array of the signal
    filter_method: callable, the method to be used for filtering
    neighbourhood_size: float, the 'diameter' of the neighbourhood to be removed around points
    padding_size=None, the size of the padding i.e. the extension of the signal. if none, padding_size = len(y) and therefore the signal is copied both ways
    points=(0, 1), the points around which the neighbourhood is to be removed
    force_zero:bool=True, whether after removing a certain neighbourhood, we force y(points) to go through 0
    filter_parameters=(, ),	 filter parameters that the filtering script may need

    1. Extends signal by padding_size
    2. Filters the extended signal
    3. Removes the points on a neighbourhood of the boundaries (0, 1)
     3a. If force_values=True, force the signal to go through the forced values (initially 0)
    4. Interpolates over the void area
    5. Returns the signal only in the interval (0, 1)

    """
    # First extend the signal
    x_extended, y_extended = extend_signal(x, y, padding_size=padding_size)

    # Second filter it
    y_extended_filtered = filter_rho(y_extended)
    y_extended_filtered = filter_rho(y_extended_filtered)

    # Third and fourth remove the boundaries and interpolate
    x_extended_removed, y_extended_filtered_removed = remove_and_interpolate(
        x_extended,
        y_extended_filtered,
        points=points,
        neighbourhood_size=neighbourhood_size,
        force_zero=force_zero,
    )

    # Fifth return the signal only in the interval 0 to 1
    x_0_to_1, y_0_to_1 = return_to_0_1(x_extended_removed, y_extended_filtered_removed)

    return x_0_to_1, y_0_to_1

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

def generate_solution_family_guess(lambda_value, max_N):
    lambda_value = lambda_value # From the way the guess functions work
    eigenvalue_guess_array = [
            eigenvalue_guess(n) for n in range(-max_N, max_N+1)
            ]
    eigenstate_guess_array = [
            dirichlet_eigenstate(z, omega_n) 
            for omega_n in eigenvalue_guess_array
            ]
    eigenstate_gradient_guess_array = [
            dirichlet_eigenstate_gradient(z, omega_n) 
            for omega_n in eigenvalue_guess_array
            ]

    solution_family = {
            "eigenvalues":eigenvalue_guess_array,
            "eigenstates":eigenstate_guess_array,
            "eigenstate_gradients":eigenstate_gradient_guess_array,
            }

    return solution_family

def save_solutions(
        converged_dictionary,
        directory="",
        sig_digs=3,
        ):
    """
    Saves the calculated quantities to {directory}/{boundary_conditions}/{quantity}/a_{a}_mass_{mass}_lambda_value_{lambda_value}.csv
    Parameters:
        solution_family: A dictionary with keys() = ["eigenvalues", "eigenstates", "eigenstate_gradients"]
        directory: In case a further directory should be considered, e.g. if Ambjorn technique is used
    Returns None
    """
    lambda_string = float_to_str(lambda_value, sig_digs=sig_digs)
    a_string = float_to_str(a, sig_digs=sig_digs)
    m_string = float_to_str(m, sig_digs=sig_digs)

    file_id = f"mass_{m_string}_a_{a_string}_lambda_{lambda_string}.txt"

    if directory != "":
        directory = "/" + directory

    root_directory = f"saved_solutions{directory}/{bcs}"
    print(f"Saving results under {root_directory}/.../{file_id}...")
    
    try:
        for key, value in converged_dictionary.items():
            savetxt(
                    f"{root_directory}/{key}/{file_id}",
                value,
                delimiter=",",
            )
    except FileNotFoundError as e:
        print(e)
        create = input(f"\nCreate directory {directory[1:]}?[y/n]... ")
        if create == "y":
            for key in converged_dictionary.keys():
                os.makedirs(root_directory + "/" + key)
            save_solutions(converged_dictionary, directory=directory, sig_digs=sig_digs)
        elif create == "n": 
            rename = input(f"If {directory} was a typo, enter the correct name...")
            if rename != "":
                save_solutions(converged_dictionary, directory=directory, sig_digs=sig_digs)



def read_solutions_from_file(
        directory="",
        sig_digs=3,
        ):
    m_string = float_to_str(m, sig_digs)
    lambda_string = float_to_str(lambda_value, sig_digs)
    a_string = float_to_str(a, sig_digs)

    file_id = f"mass_{m_string}_a_{a_string}_lambda_{lambda_string}.txt"

    if directory != "": 
        directory = "/" + directory

    # Initialize the solution family dict
    solution_family = {
            "eigenvalues":0,
            "eigenstates":0,
            "eigenstate_gradients":0,
            "rho":0,
            "A0_induced":0,
            }
    
    for key in solution_family.keys():
        solution_family[key]  = np.genfromtxt(
            f"saved_solutions{directory}/{bcs}/{key}/{file_id}",
            dtype=float,
            delimiter="\n",
        )

    return solution_family

def single_system_iteration():
    global rho, solution_family, A0_induced, A0
    solution_family = calculate_eigenstates(
            max_N,
            max_nodes=max_nodes,
            tol=tol,
            )

    solution_family = normalize_eigenstate_family(solution_family)

    rho = calculate_rho(solution_family)
    if smoothing:
        # rho = filter_rho(rho)
        rho += e/np.pi * A0(z)
        _, rho = extend_and_filter(z, rho, )

    A0_induced = calculate_A0_induced(rho)
    A0 = calculate_A0_total(A0_induced)

def update_system_iterate(n_iterations, plot, r=1):
    for n in count(0):
        if isinstance(n_iterations, int):
            if n >= n_iterations: break
        elif n_iterations is None:
            if r < iterate_tol: break
        else: 
            raise ValueError("n_iterations must be either int or None")

        single_system_iteration()

        if plot and not i%3:
            plot_rho_A0_induced(rho, A0_induced(z), '--', color=color, alpha=0.5)
        if not n: # For the first iteration
            print('n=0, max(A0_induced) =', max(A0_induced(z)))
        
        if n_iterations is None:
            try:
                r = abs(max(rho) - max_prev_rho) / abs(max_prev_rho)
                max_prev_rho = max(rho)
            except NameError:
                max_prev_rho = max(rho)
    return n

if __name__ == "__main__":
    # When wanting to calculate with n_points = 100, max_N = 12
    # it is the return_to_0_1 what gives troubles. I don't see why
     
    plot = True
    save = False
    read = False

    if plot:
        fig, (ax_rho, ax_A0_induced, ax_eigenvalues) = plt.subplots(3)
        eigenvalues_array = []
        lambda_array = []

    # Number of modes
    n_points = 88
    # x-axis
    z = np.linspace(0, 1, n_points)

    e = 1
    a = 5
    m = 0
    max_N = 10 # Cutoff for the mode solutions
               # Note that there will be 2 * max_N - 1 solutions
    lambda_min = 7 # Min strength of the external electric field
    lambda_max = 25 # Max strength of the external electric field
    lambda_step= 0.25

    tol = 1e-10 # Tolerance for solve_bvp
    iterate_tol = 1e-7 # Tolerance for convergence of selfconsistent solutions

    max_nodes = 5e23
    smoothing = True
    ambjorn = False
    bcs = "dirichlet"
    n_iterations = 2 # Number of iterations for each value of lambda
    directory = "n_iterations_2_hadamard_maxN_10"

    print(f"max_N = {max_N}")
    print(f"smoothing = {smoothing}")
    print(f"ambjorn = {ambjorn}")
    print(f"n_iterations = {n_iterations}")
    if save:
        print(f"Saving to directory {directory}")
    else:
        print("Not saving the solutions")

    # Need to initialize the A0_induced
    A0_induced = CubicSpline(z, np.zeros_like(z))
    color_cycle = cycle(["#FE5F55","#3423A6","#710627","#7180B9","#163B20"])

    lambda_value = lambda_min

    if read:
        try: 
            raise Exception("Solutions can't be read yet, eigenstates aren't read correctly")
            solution_family = read_solutions_from_file(
                    directory=directory,
                    )
            A0_induced = solution_family["A0_induced"]
        except FileNotFoundError:
            solution_family = generate_solution_family_guess(lambda_min, max_N)
            A0_induced = CubicSpline(z, np.zeros_like(z))
            pass
    else:
        solution_family = generate_solution_family_guess(lambda_min, max_N)
        A0_induced = CubicSpline(z, np.zeros_like(z))

    for i, lambda_value in enumerate(np.arange(lambda_min, lambda_max, lambda_step)):
        r = 1 # In case n_iterations is None
        print(15*"#")
        print('λ= ', lambda_value)
        color = next(color_cycle) # To distinguish between different lambda values
        try:
            n = update_system_iterate(n_iterations, plot) # The number of iterations needed to converge
            if plot:
                eigenvalues_array.append(solution_family["eigenvalues"])
                lambda_array.append(lambda_value)
                if not i%3:
                    plot_rho_A0_induced(rho, A0_induced(z), '', color=color, alpha=1)

            print(f'n={n}, max(A0_induced) =', max(A0_induced(z)))
            if save:
                solution_family["rho"] = rho
                solution_family["A0_induced"] = A0_induced(z)
                save_solutions(
                    solution_family,
                    directory=directory,
                    )

            A0 = calculate_A0_total(lambda z: (lambda_value + lambda_step)/lambda_value * A0_induced(z))

        except Exception as exception:
            ex = exception
            break
        except KeyboardInterrupt as exception:
            ex = exception
            break

    if plot:
        ax_eigenvalues.plot(lambda_array, eigenvalues_array, 'b')
        plt.show()
    
    raise ex
