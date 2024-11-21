import numpy as np
import os
from itertools import cycle
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from pathlib import Path
from mpmath import quad

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

class Vacuum_Polarization:
    def __init__(self,
            max_N=1,
            m=0,
            a=1,
            e=1,
            lambda_min=10,
            lambda_max=25,
            lambda_step=2,
            relax_parameter=0.2,
            max_n_iterations=200,
            walkback_factor=0.4,
            lambda_step_min=1e-6,
            n_iterations=None,
            iterate_tol=1e-3,
            smoothing=True,
            ambjorn=True,
            plot = True, 
            plot_for_each_lambda = True, 
            save = True, 
            save_plots = True, 
            show_plots = False, 
            read = False, 
            directory = "",
            bcs = "dirichlet",
            ):

        self.e = e
        self.a = a
        self.max_N = max_N
        self.m = m
        self.eigenvalues = [ 
            self.eigenvalue_guess(n) for n in range(-max_N, max_N+1) if n!=0
            ]

        self.color_cycle = cycle(["#FE5F55","#3423A6","#710627","#7180B9","#163B20"])
         
        self.plot = plot
        self.plot_for_each_lambda = plot_for_each_lambda
        self.save = save
        self.save_plots = save_plots
        self.show_plots = show_plots
        self.read = read

        self.plotting_setup()

        self.relax_parameter = relax_parameter # the number c in A_{k+1} = (1-c) A_k + c ( -𝜆 (z-1/2) - ∫∫𝜌)

        if max_N<5:
            self.n_points = 200
            self.ambjorn = True # This can't be changed
            self.smoothing = False
        else:
            self.n_points = 8 * ( max_N + 1 )
            self.ambjorn = ambjorn # This can be changed depending on the results we want
            self.smoothing = smoothing # This can be changed depending on the results we want


        self.z = np.linspace(0, 1, self.n_points)
        
        self.lambda_min = lambda_min
        self.lambda_max = lambda_max
        self.lambda_step = lambda_step
        self.lambda_value = lambda_min

        self.A0_induced = lambda z: 0
        self.A0_induced_history = [[self.A0_induced]]
        self.A0 = lambda z: - self.lambda_value * (z - 1/2) 
        self.A0_history = [[self.A0]]

        self.n_iterations = n_iterations # In case I want it to just iterate a small finite amount of times
                    # if n_iterations = None, "back-react" until convergence
        self.iterate_tol = iterate_tol # Tolerance for the convergence of the backreaction procedure
        self.max_n_iterations = max_n_iterations # Avoid the calculations from getting in a loop
        self.walkback_factor = walkback_factor # Avoid the calculations from getting in a loop
        self.lambda_step_min = lambda_step_min # Avoid the calculations from getting in a loop

        self.directory = directory
        self.bcs = bcs

        if save_plots:

            directory_exists = Path("figures/"+directory+"/"+bcs).is_dir()

            if directory_exists:
                overwrite = input("You are about to overwrite an existing data folder. Continue? [y/n]")
                if overwrite=="y":
                    shutil.rmtree("figures/"+directory+"/"+bcs)
                elif overwrite=="n":
                    print("Then change the code")
                    exit()

            Path("figures/"+directory+"/"+bcs).mkdir(parents=True, exist_ok=True)


    def plotting_setup(self):
        global plt
        import matplotlib.pyplot as plt
        if self.plot:
            self.fig1 = plt.figure(r"Vacuum polarization")
            self.fig2 = plt.figure(r"A_0 induced")
            self.fig3 = plt.figure("Mode energy evolution")
            self.fig4 = plt.figure("Values at intermediate steps")

            self.ax_rho = self.fig1.subplots(1)
            self.ax_A0_induced = self.fig2.subplots(1)
            self.ax_eigenvalues = self.fig3.subplots(1)
            self.ax_relax, self.ax_intermediate_omegas = self.fig4.subplots(2)

            self.eigenvalues_array = []
            self.lambda_array = []

    # Some util functions
    def sign(self, x):
        # Returns sign of x
        return 2 * (x>0) - 1 

    def float_in_array(self, x, array, tol=1e-3): # This is not used
        # Checks if x is in given array up to some tolerance
        for element in array:
            if abs(x - element) < tol: 
                return True
        return False

    def float_to_str(self, value, sig_digs=3):
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

    def str_to_float(self, value:str):
        return float(
            value.replace("_", ".")
            )

    def root_mean_square(self, y):
        return np.sqrt(np.mean(y**2))

    def find_root_bisection(self, fun, x1, x2, tol=1e-4, maxiter=500):
        """Uses the bisection method to find the root of a function. Works recursively
        Parameters:
            fun: the callable of which to find the roots
            x1: a guess for the root
            x2: a second guess for the root
            tol=1e-4: the precision with which the roots are searched
            maxiter=500: the maximum amount of iterations before saying no solutions were found
        Returns:
            c: the solution
        """
        if abs(fun(x1)) < tol:
            return x1
        elif abs(fun(x2)) < tol:
            return x2

        if maxiter==0:
            raise RuntimeError(f"maxiter was reached")
        if fun(x1)*fun(x2) > 0:
            raise RuntimeError(f"No root was found in the interval ({x1}, {x2})")

        c = (x1 + x2) / 2
        if fun(x1)*fun(c) < 0: # Means the root was at the left of c
            return self.find_root_bisection(fun, x1, c, tol=tol, maxiter=maxiter-1)
        else:
            return self.find_root_bisection(fun, c, x2, tol=tol, maxiter=maxiter-1)

    def bisection_method_upper_bound(self, previous_array):
        """
        Given an array of guesses for the roots of a function,
        returns an array of upper bound guesses
        """
        if len(previous_array) == 2: return [1.9 * element for element in previous_array] 
        
        positive_solutions = [element for element in previous_array if element > 0]
        guess = [0] * len(positive_solutions)

        for index, element in enumerate(positive_solutions[:-1]):
            guess[index]= (element + positive_solutions[index+1])/2

        # Some obscure math was used  to get to this value.
        guess[-1] = (3 * positive_solutions[-1] - positive_solutions[-2])/2

        return [-element for element in guess[::-1]] + guess

    def bisection_method_lower_bound(self, previous_array):
        if len(previous_array) == 2: return [0., 0.]
        
        positive_solutions = [element for element in previous_array if element > 0]
        guess = [0] * len(positive_solutions)

        for index, element in enumerate(positive_solutions):
            if index==0:
                continue

            guess[index]= (element + positive_solutions[index-1])/2

        return [-element for element in guess[::-1]] + guess

    def save_solutions(
            self,
            sig_digs=3,
            ):
        """
        Saves the calculated quantities to {directory}/{boundary_conditions}/{quantity}/mass_{mass}_a_{a}_lambda_value_{lambda_value}.csv
        Parameters:
            solution_family: A dictionary with keys() = ["eigenvalues", "eigenstates", "eigenstate_gradients"]
            directory: In case a further directory should be considered, e.g. if Ambjorn technique is used
        Returns None
        """

        solution_family = {
                "eigenvalues":self.eigenvalues,
                "eigenstates":self.eigenstates,
                "A0_induced":self.A0_induced(self.z),
                "rho":self.rho,
                }

        lambda_string = self.float_to_str(self.lambda_value, sig_digs=sig_digs)
        a_string = self.float_to_str(self.a, sig_digs=sig_digs)
        m_string = self.float_to_str(self.m, sig_digs=sig_digs)

        file_id = f"mass_{m_string}_a_{a_string}_lambda_{lambda_string}.txt"

        if self.directory != "":
            directory = "/" + self.directory

        root_directory = f"saved_solutions{directory}/{self.bcs}"
        print(f"Saving results under {root_directory}/.../{file_id}...")
        
        try:
            for key, value in solution_family.items():
                np.savetxt(
                        f"{root_directory}/{key}/{file_id}",
                    value,
                    delimiter=",",
                )
        except FileNotFoundError as e:
            print(e)
            create = input(f"\nCreate directory {directory[1:]}?[y/n]... ")
            if create == "y":
                for key in solution_family.keys():
                    os.makedirs(root_directory + "/" + key)
                self.save_solutions()
            elif create == "n": 
                rename = input(f"If {directory} was a typo, enter the correct name...")
                if rename != "":
                    self.save_solutions()

    def eigenvalue_guess(self, n):
        return self.sign(n) * np.sqrt( (n * np.pi)**2 + self.m ** 2 )

    def Klein_Gordon_equation( self, z_array, y, omega):
        # The differential equation

        background_field = self.A0(z_array)

        klein_gordon = np.array(
            (y[1], -((omega - background_field) ** 2 + self.m ** 2) * y[0])
        )

        return klein_gordon


    def calculate_eigenstates(self):
        """
        Calculates the solutions to the klein_gordon equation by parametrizing the solution family by a parameter, and then looking for the parameter that verifies the boundary conditions
        """

        if self.bcs == "dirichlet":
            initial_values = (0, 1)
            bcs_index = 0 
        elif self.bcs == "neumann":
            initial_values = (1, 0)
            bcs_index = 1 

        eigenvalue_lower_bound = self.bisection_method_lower_bound(self.eigenvalues)
        eigenvalue_upper_bound = self.bisection_method_upper_bound(self.eigenvalues)

        parametrized_ODE = lambda omega: solve_ivp(
                lambda z, y: self.Klein_Gordon_equation(z, y, omega), (0,1), initial_values, dense_output=True
                ) # The first 0 1 is the range, the second one are the initial values

        # import matplotlib.pyplot as plt


        # omega_array = np.arange(0, 20, 0.1)
        # plt.plot(omega_array, [parametrized_ODE(omega).sol(1)[1] for omega in omega_array], )
        # plt.show()

        # exit()
        eigenvalues = [
                self.find_root_bisection(
                    lambda omega: parametrized_ODE(omega).sol(1)[bcs_index], *omega_upper_lower
                    ) 
                for omega_upper_lower in zip(
                    eigenvalue_lower_bound,
                    eigenvalue_upper_bound
                    )
                ]

        # the eigenvalues should be antisymmetric i.e. omega_n = -omega_{-n}
        self.eigenvalues = [ (i-j) / 2 for i, j in zip(eigenvalues, eigenvalues[::-1]) ]

        self.eigenstates = [ parametrized_ODE(omega).sol(self.z)[bcs_index] for omega in eigenvalues ]

    def normalize_eigenstates(self):
        """
        Normalizes the eigenstates 
        """
        # Warning, math

        for n, (eigenvalue, eigenstate) in enumerate(
                zip(
                    self.eigenvalues,
                    self.eigenstates,
                )
                ):

            eigenstate = CubicSpline(self.z, eigenstate)

            def rho_n_without_normalizing(z):
                # Normalizing wrt the symplectic norm
                # the solutions need not be real.
                return (eigenvalue - self.A0(z)) * abs(eigenstate(z))**2

            # Calculate the norm
            norm_squared = abs(float(quad(rho_n_without_normalizing, [0, 1])))

            norm = self.sign(eigenvalue)*np.sqrt(norm_squared)

            self.eigenstates[n] /= norm

    def calculate_rho(self, ):
        """
        Calculates the vacuum polarization AFTER normalizing the eigenstates, but before filtering
        """
    
        # The total charge density. 
        # Initiate as np.array to add them
        rho = np.zeros_like(self.z) 

        for n, (eigenvalue, eigenstate) in enumerate(
                zip(
                    self.eigenvalues,
                    self.eigenstates,
                    )
                ):
            # The charge density associated to the nth mode
            # eigenstate is np array. A0 is a callable. 
            # A0(z) is of the same shape as eigenstate
            rho_n = (eigenvalue - self.A0(self.z)) * abs(eigenstate) ** 2
            rho += 1/2*rho_n
    
        return rho

    def calculate_relax_parameter(self):
        """
        In the update law A_{k+1} = (1-c) A_k + c ( - λ (z - 1/2) - ∫ ∫ ρ ), 
        calculate c so that the correction term c( - λ ... ) is exactly self.correction_parameter
        """
        c = (
                self.correction_parameter * (- self.lambda_value / 2 + self.a * self.A0_induced(1) )  
                / ( (1-self.correction_parameter) * self.constant_lambda_A0_list[-1](1) 
                    + self.correction_parameter * ( - self.lambda_value / 2 + self.a * self.A0_induced(1))
                    )
                )
        return c

    def calculate_A0_induced(self, rho):
        """
        Calculates the induced A0 as the solution to the Poisson equation.
        It is just integrating the vacuum polarization twice
        """

        rho_interpolated = CubicSpline(self.z, self.rho)
        A0_induced_shifted = rho_interpolated.antiderivative(2)
        offset = A0_induced_shifted(1/2)
        A0_induced = lambda z: -(A0_induced_shifted(z) - offset)

        if A0_induced(1) > 10:
            print(f"n={self.n}, A0_induced(1)={self.A0_induced(1)}")

        return A0_induced

    def plot_rho_A0_induced(self, rho, A0_induced, fmt, color, alpha,):
        self.ax_rho.plot(self.z, rho, fmt, color=self.color, alpha=alpha)
        self.ax_A0_induced.plot(self.z, A0_induced, fmt, color=self.color, alpha=alpha)

    def single_system_iteration(self):

        self.calculate_eigenstates()
        self.normalize_eigenstates()

        rho_new = self.calculate_rho()
        if self.smoothing and not self.ambjorn:
            # rho = filter_rho(rho)
            rho_new += e**2/np.pi * self.A0(self.z_array)
            _, rho_new = extend_and_filter(self.z_array, self.rho)

        if self.n>70 and not self.n%25 and False:
            print("Averaging..")
            self.rho = (rho_new + self.rho)/2 
        else:
            self.rho = rho_new

        self.A0_induced = self.calculate_A0_induced(self.rho)

    def update_system_iterate(self):
        r=1

        if self.save_plots and self.plot_for_each_lambda:

            self.fig4.suptitle(r"$\lambda$={}".format(round(self.lambda_value, 4)))

            self.ax_relax.cla()
            self.ax_intermediate_omegas.cla()

            self.ax_intermediate_omegas.set_xlabel(r"n")
            self.ax_relax.set_ylabel(r"$A_0(1) +$ {}".format(round(self.lambda_value/2, 5)))
            self.ax_intermediate_omegas.set_ylabel(r"$\omega$")

        self.n=0

        self.constant_lambda_A0_list = [ self.A0_history[-1][-1] ]
        self.constant_lambda_A0_induced_list = [ self.A0_induced_history[-1][-1] ]

        while self.n < self.max_n_iterations:

            if self.relax_parameter is None:
                c = self.calculate_relax_parameter()
            else:
                c = self.relax_parameter

            self.A0_induced = self.constant_lambda_A0_induced_list[-1]
            self.A0 = lambda z:  (1-c) * self.constant_lambda_A0_list[-1](z) + c * (- self.lambda_value * (z-1/2) + self.a * self.A0_induced(z))

            if isinstance(self.n_iterations, int):
                if self.n >= self.n_iterations: break
            elif self.n_iterations is None:
                if r < self.iterate_tol: break
            else: 
                raise Exception("n_iterations must be either int or None")

            self.single_system_iteration()

            if self.plot:
                self.plot_rho_A0_induced(self.rho, self.A0_induced(self.z), '-', color=self.color, alpha=0.2)

            if self.plot_for_each_lambda and self.save_plots:
                self.ax_relax.plot(self.n, self.A0(1) + round(self.lambda_value/2, 2), 'x', color=self.color)
                self.ax_intermediate_omegas.plot(self.n, self.eigenvalues[self.max_N], 'x', color=self.color)

            if self.n == 0: # For the first iteration
                print('n=0, A0_induced(1) =', self.A0_induced(1))
            
            if self.n_iterations is None:
                try:
                    r = abs(max(self.rho) - max_prev_rho) / abs(max_prev_rho)
                    max_prev_rho = max(self.rho)
                except NameError:
                    max_prev_rho = max(self.rho)

            self.A0 = CubicSpline(self.z, self.A0(self.z))

            self.constant_lambda_A0_list.append(self.A0)
            self.constant_lambda_A0_induced_list.append(self.A0_induced)
            self.n += 1

        else:
            if self.save_plots:
                self.fig4.savefig("figures/"+self.directory+"/"+self.bcs+"/intermediate_steps_lambda_value_"+self.float_to_str(self.lambda_value)+".png")
            raise RuntimeError("max_n_iterations reached")

        if self.save_plots: 
            self.fig4.savefig("figures/"+self.directory+"/"+self.bcs+"/intermediate_steps_lambda_value_"+self.float_to_str(self.lambda_value)+".png")

    def full_script(self):

        while self.lambda_value < self.lambda_max:
            print(self.eigenvalues[self.max_N-2:self.max_N+2])

            if self.lambda_step < self.lambda_step_min:
                break

            print(15*"#")
            print('λ =', self.lambda_value)
            self.color = next(self.color_cycle) # To distinguish between different lambda values

            try:
                self.update_system_iterate() # The number of iterations needed to converge
                if self.plot:
                    self.eigenvalues_array.append(self.eigenvalues)
                    self.lambda_array.append(self.lambda_value)
                    self.plot_rho_A0_induced(self.rho, self.A0_induced(self.z), '', color=self.color, alpha=1)

                print(f'Converged in n={self.n} iterations. A0_induced(1) =', self.A0_induced(1))
                if self.save:
                    solution_family = {"eigenvalues":self.eigenvalues}
                    solution_family["rho"] = self.rho
                    solution_family["A0_induced"] = self.A0_induced(self.z)
                    self.save_solutions()

            except RuntimeError as exception:
                self.lambda_value -= self.lambda_step 
                self.lambda_step *= self.walkback_factor
                self.relax_parameter *= self.walkback_factor
                self.lambda_value += self.lambda_step

                print(exception) # To know if it didn't converge or if no solution was found
                print(f'Error found at n={self.n}, A0_induced(1) =', self.A0_induced(1))

                continue
            except Exception as exception:
                ex = exception
                break
            except KeyboardInterrupt as exception:
                ex = exception
                break

            self.A0_history.append(self.constant_lambda_A0_list)
            self.A0_induced_history.append(self.constant_lambda_A0_induced_list)


            self.lambda_value += self.lambda_step

        
        if self.save_plots:
            self.ax_eigenvalues.plot(self.lambda_array, self.eigenvalues_array, 'b')
            self.fig1.savefig("figures/"+self.directory+"/"+self.bcs+"/vacuum_polarization.png")
            self.fig2.savefig("figures/"+self.directory+"/"+self.bcs+"/A0_induced.png")
            self.fig3.savefig("figures/"+self.directory+"/"+self.bcs+"/eigenvalues.png")
        if self.show_plots:
            self.ax_eigenvalues.plot(self.lambda_array, self.eigenvalues_array, 'b')
            plt.show()

        raise ex
