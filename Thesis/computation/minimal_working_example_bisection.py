from pathlib import Path
import shutil
from mpmath import quad
from itertools import count, product
import numpy as np
from numpy import savetxt
from scipy.integrate import solve_ivp
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
    return - lambda_value * (z- 1/2) + a * A0_induced(z)

def eigenvalue_guess(n):
    return n * np.pi

def Klein_Gordon_equation(
    z_array,
    y,
    omega
    ):
    # The differential equation

    background_field = A0(z_array)

    klein_gordon = np.array(
        (y[1], -((omega - background_field) ** 2 + m ** 2) * y[0])
    )

    return klein_gordon

def find_root_bisection(fun, x1, x2, tol=1e-4, maxiter=500):
    if abs(fun(x1)) < tol:
        return x1
    elif abs(fun(x2)) < tol:
        return x2

    if maxiter==0:
        raise RecursionError(f"maxiter was reached")
    if fun(x1)*fun(x2) > 0:
        raise ValueError(f"No root was found in the interval ({x1}, {x2})")

    c = (x1 + x2) / 2
    if fun(x1)*fun(c) < 0: # Means the root was at the left of c
        return find_root_bisection(fun, x1, c, tol=tol, maxiter=maxiter-1)
    else:
        return find_root_bisection(fun, c, x2, tol=tol, maxiter=maxiter-1)

# Since bisection method needs two starting points, here I generate
# lower and upper bounds of the eigenvalues. They should always be in these 
# intervals
def bisection_method_upper_bound(previous_array):
    if len(previous_array) == 2: return [-3.5, 3.5] 
    
    positive_solutions = [element for element in previous_array if element > 0]
    guess = [0] * len(positive_solutions)

    for index, element in enumerate(positive_solutions[:-1]):
        guess[index]= (element + positive_solutions[index+1])/2

    # Some obscure math was used  to get to this value.
    guess[-1] = (3 * positive_solutions[-1] - positive_solutions[-2])/2

    return [-element for element in guess[::-1]] + guess

def bisection_method_lower_bound(previous_array):
    if len(previous_array) == 2: return [0., 0.]
    
    positive_solutions = [element for element in previous_array if element > 0]
    guess = [0] * len(positive_solutions)

    for index, element in enumerate(positive_solutions):
        if index==0:
            continue

        guess[index]= (element + positive_solutions[index-1])/2

    return [-element for element in guess[::-1]] + guess

def calculate_eigenstates():

    eigenvalues = []
    eigenstates = []

    eigenvalue_guess_array = solution_family["eigenvalues"]

    eigenvalue_lower_bound = bisection_method_lower_bound(eigenvalue_guess_array)
    eigenvalue_upper_bound = bisection_method_upper_bound(eigenvalue_guess_array)

    parametrized_ODE = lambda omega: solve_ivp(
            lambda z, y: Klein_Gordon_equation(z, y, omega), (0,1), (0,1), dense_output=True
            ) # The first 0 1 is the range, the second one are the initial values

    eigenvalues = [
            find_root_bisection(
                lambda omega: parametrized_ODE(omega).sol(1)[0], *omega_upper_lower
                ) 
            for omega_upper_lower in zip(
                eigenvalue_lower_bound,
                eigenvalue_upper_bound
                )
            ]

    eigenvalues = [(i-j) / 2 for i, j in zip(eigenvalues, eigenvalues[::-1])] # the eigenvalues should be antisymmetric i.e. omega_n = -omega_{-n}

    eigenstates = [ parametrized_ODE(omega).sol(z)[0] for omega in eigenvalues ]

    return {
            "eigenvalues": eigenvalues,
            "eigenstates": eigenstates
            }

def normalize_eigenstate_family(solution_family:dict):
    # Due to the boundary conditions, the solutions are not normalized
    # This assumes that the A0 is somewhere defined.

    # Had to initiate the dict again, I am not sure why
    # If I don't do it, it doesn't really update.
    eigenvalues = solution_family["eigenvalues"]
    eigenstates = []
    # Warning, math
    for n, (eigenvalue, eigenstate) in enumerate(
            zip(
            solution_family["eigenvalues"],
            solution_family["eigenstates"],
            )
            ):

        eigenstate = CubicSpline(z, eigenstate)

        def rho_n_without_normalizing(z):
            # Normalizing wrt the symplectic norm
            # the solutions need not be real.
            return (eigenvalue - A0(z)) * abs(eigenstate(z))**2

        # Calculate the norm
        norm_squared = abs(float(quad(rho_n_without_normalizing, [0, 1])))

        norm = sign(eigenvalue)*np.sqrt(norm_squared)

        eigenstates.append(eigenstate(z)/norm)

    return {
            "eigenvalues":eigenvalues,
            "eigenstates":eigenstates,
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

    if A0_induced(1) > 10:
        print(f"n={n}, A0_induced(1)={A0_induced(1)}")

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

def generate_solution_family_guess(lambda_value, max_N):
    lambda_value = lambda_value # From the way the guess functions work
    eigenvalue_guess_array = [
            eigenvalue_guess(n) for n in range(-max_N, max_N+1) if n!=0
            ]

    solution_family = {
            "eigenvalues":eigenvalue_guess_array,
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
            )

    solution_family = normalize_eigenstate_family(solution_family)

    rho_new = calculate_rho(solution_family)
    if smoothing and not ambjorn:
        # rho = filter_rho(rho)
        rho_new += e**2/np.pi * A0(z)
        _, rho_new = extend_and_filter(z, rho_new)

    if n>70 and not n%25:
        rho = (rho_new + rho)/2 
    else:
        rho = rho_new

    A0_induced = calculate_A0_induced(rho)

def update_system_iterate(n_iterations, plot, r=1):
    if plot_for_each_lambda:

        fig4.suptitle(r"$\lambda$={}".format(round(lambda_value, 4)))

        ax_relax.cla()
        ax_intermediate_omegas.cla()

        ax_intermediate_omegas.set_xlabel(r"n")
        ax_relax.set_ylabel(r"$A_0(1) +$ {}".format(round(lambda_value/2, 5)))
        ax_intermediate_omegas.set_ylabel(r"$\omega$")

    global A0, n
    n=0

    while n < max_n_iterations:

        if relax_parameter is None:
            c = (
                    correction_parameter * (- lambda_value / 2 + a * A0_induced(1) )  
                    / ( (1-correction_parameter) * constant_lambda_A0_list[-1](1) 
                        + correction_parameter * ( - lambda_value / 2 + a * A0_induced(1))
                        )
                    )
        else:
            c = relax_parameter

        A0 = lambda z:  (1-c) * constant_lambda_A0_list[-1](z) + c * (- lambda_value * (z-1/2) + a * A0_induced(z))

        if isinstance(n_iterations, int):
            if n >= n_iterations: break
        elif n_iterations is None:
            if r < iterate_tol: break
        else: 
            raise Exception("n_iterations must be either int or None")

        single_system_iteration()

        if plot:
            plot_rho_A0_induced(rho, A0_induced(z), '-', color=color, alpha=0.2)

        if plot_for_each_lambda:
            ax_relax.plot(n, A0(1) + round(lambda_value/2, 2), 'x', color=color)
            ax_intermediate_omegas.plot(n, solution_family["eigenvalues"][max_N], 'x', color=color)

        if not n: # For the first iteration
            print('n=0, A0_induced(1) =', A0_induced(1))
        
        if n_iterations is None:
            try:
                r = abs(max(rho) - max_prev_rho) / abs(max_prev_rho)
                max_prev_rho = max(rho)
            except NameError:
                max_prev_rho = max(rho)


        A0 = CubicSpline(z, A0(z))

        constant_lambda_A0_list.append(A0)
        n += 1

    else:
        if save_plots:
            fig4.savefig("figures/"+directory+"/"+bcs+"/intermediate_steps_lambda_value_"+float_to_str(lambda_value)+".png")
        raise ValueError("max_n_iterations reached")

    if save_plots:
        fig4.savefig("figures/"+directory+"/"+bcs+"/intermediate_steps_lambda_value_"+float_to_str(lambda_value)+".png")

    return n

if __name__ == "__main__":
    # When wanting to calculate with n_points = 100, max_N = 12
    # it is the return_to_0_1 what gives troubles. I don't see why
    color_cycle = cycle(["#FE5F55","#3423A6","#710627","#7180B9","#163B20"])
     
    plot = True
    plot_for_each_lambda = True
    save = True
    save_plots = True
    show_plots = False
    read = False

    if save_plots or plot:
        import matplotlib.pyplot as plt
        fig1 = plt.figure(r"Vacuum polarization")
        ax_rho = fig1.subplots(1)
        fig2 = plt.figure(r"A_0 induced")
        ax_A0_induced = fig2.subplots(1)
        fig3 = plt.figure("Mode energy evolution")
        ax_eigenvalues = fig3.subplots(1)
        fig4 = plt.figure("Values at intermediate steps")
        ax_relax, ax_intermediate_omegas = fig4.subplots(2)

        eigenvalues_array = []
        lambda_array = []

    # Number of modes
    # x-axis

    e = 1
    a = 1
    m = 0
    max_N = 12 # Cutoff for the mode solutions
               # Note that there will be 2 * max_N - 1 solutions
    if max_N<5:
        n_points = 200
        ambjorn = True # This can't be changed
    else:
        n_points = 8 * ( max_N + 1 )
        ambjorn = False # This can be changed depending on the results we want

    z = np.linspace(0, 1, n_points)

    lambda_min = 10 # Min strength of the external electric field
    lambda_max = 25 # Max strength of the external electric field
    lambda_step= 1
    relax_parameter = 0.2
    correction_parameter = 0.95 # So that the correction is only a 10% of the existing potential

    walkback_factor = 0.3
    lambda_step_min = 1e-6

    iterate_tol = 1e-3 # Tolerance for convergence of selfconsistent solutions

    smoothing = True
    bcs = "dirichlet"
    n_iterations = None # Number of iterations for each value of lambda
    max_n_iterations = 200 # break convergence iterations if 
    directory = "hadamard_averaging"

    print(f"max_N = {max_N}")
    print(f"smoothing = {smoothing}")
    print(f"ambjorn = {ambjorn}")
    print(f"n_iterations = {n_iterations}")
    print("plot =", plot)
    print("plot_for_each_lambda =", plot)
    print("save =",save)
    print("save_plots =",save_plots)
    print("show_plots =",show_plots)
    print("read =",read)
    if save:
        print(f"Saving to directory {directory}")
    else:
        print("Not saving the solutions")

    if save_plots:

        directory_exists = Path("figures/"+directory+"/"+bcs).is_dir()

        if directory_exists:
            overwrite = input("You are about to overwrite an existing data folder. Continue? [y/n]")
            if overwrite=="y":
                shutil.rmtree("figures/"+directory+"/"+bcs)
            elif overwrite=="n":
                save_plots = input("Do you want to save plots? [y/n]")
                if save_plots == "n":
                    save_plots=False
                elif save_plots == "y":
                    directory = input("What should the name of the directory be instead?")

        Path("figures/"+directory+"/"+bcs).mkdir(parents=True, exist_ok=True)

    # Need to initialize the A0_induced
    A0_induced = lambda z: 0

    lambda_value = lambda_min

    A0_list = [[ lambda z: -lambda_value * (z-1/2) ]]

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


    # Main calculation
    i=0
    while lambda_value < lambda_max:
        print(solution_family["eigenvalues"][max_N-2:max_N+2])

        if lambda_step < lambda_step_min:
            break

        print(15*"#")
        print('λ= ', lambda_value)
        color = next(color_cycle) # To distinguish between different lambda values

        constant_lambda_A0_list = [A0_list[-1][-1]]

        A0_induced = lambda z: constant_lambda_A0_list[-1](z) + lambda_value * (z - 1/2)

        try:
            n = update_system_iterate(n_iterations, plot) # The number of iterations needed to converge
            if plot:
                eigenvalues_array.append(solution_family["eigenvalues"])
                lambda_array.append(lambda_value)
                plot_rho_A0_induced(rho, A0_induced(z), '', color=color, alpha=1)

            print(f'Converged in n={n} iterations. A0_induced(1) =', A0_induced(1))
            if save:
                solution_family["rho"] = rho
                solution_family["A0_induced"] = A0_induced(z)
                save_solutions(
                    solution_family,
                    directory=directory,
                    )

        except ValueError as exception:
            lambda_value -= lambda_step 
            lambda_step *= walkback_factor
            relax_parameter *= walkback_factor
            lambda_value += lambda_step

            print(exception) 
            print(f'Error found at n={n}, A0_induced(1) =', A0_induced(1))

            continue
        except Exception as exception:
            ex = exception
            break
        except KeyboardInterrupt as exception:
            ex = exception
            break

        A0_list.append(constant_lambda_A0_list)

        lambda_value += lambda_step
        # i += 1

    
    if save_plots:
        ax_eigenvalues.plot(lambda_array, eigenvalues_array, 'b')
        fig1.savefig("figures/"+directory+"/"+bcs+"/vacuum_polarization.png")
        fig2.savefig("figures/"+directory+"/"+bcs+"/A0_induced.png")
        fig3.savefig("figures/"+directory+"/"+bcs+"/eigenvalues.png")
    if show_plots:
        ax_eigenvalues.plot(lambda_array, eigenvalues_array, 'b')
        plt.show()

    
    raise ex
