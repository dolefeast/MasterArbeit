import numpy as np
import os
from itertools import cycle
from scipy.interpolate import CubicSpline
from pathlib import Path
from mpmath import quad

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
            ambjorn=False,
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

        self.color_cycle = cycle(
                [
    "#FF5733",  # Red-Orange
    "#FF8D1A",  # Golden Yellow
    "#FFD700",   # Yellow
    "#B4FF00",   # Lime Green
    "#00E5FF",   # Light Cyan
    "#1A74FF",   # Sky Blue
    "#9C1AFF",   # Vivid Purple
    "#FF007F",   # Deep Pink
    "#8000FF",   # Purple
    "#FF00FF",   # Magenta
    "#00FF00",   # Green
    "#00FFFF",   # Cyan
    "#FF6600"    # Orange
]

                )

         
        self.plot = save_plots or show_plots
        self.plot_for_each_lambda = plot_for_each_lambda
        self.save = save
        self.save_plots = save_plots
        self.show_plots = show_plots
        self.read = read

        if self.plot:
            import matplotlib.pyplot as plt
            self.plt = plt
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

        self.rho = [0] * self.n_points
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
    from Vacuum_Polarization.filter_routines import (
            extend_signal, 
            remove_neighbourhood,
            remove_and_interpolate,
            return_to_0_1,
            filter_rho,
            extend_and_filter 
            )

    from Vacuum_Polarization.small_routines import (
            sign,
            float_to_str,
            str_to_float,
            root_mean_square,
            set_config_from_dict
            )

    from Vacuum_Polarization.plotting_setup import (
            plotting_setup,
            plot_rho_A0_induced,
            plot_intermediate_steps_setup,
            plot_intermediate_steps,
            save_all_plots
            )
    from Vacuum_Polarization.bisection_routines import (
            find_root_bisection,
            bisection_method_upper_bound,
            bisection_method_lower_bound,
            )
    from Vacuum_Polarization.physics import (
            eigenvalue_guess,
            Klein_Gordon_equation
            )
    from Vacuum_Polarization.calculate_eigenstates_routines import (
            calculate_eigenstates,
            normalize_eigenstates,
            )
    from Vacuum_Polarization.save_solutions_routine import save_solutions
    from Vacuum_Polarization.induced_fields_routines import (
            calculate_rho,
            calculate_relax_parameter,
            calculate_A0_induced
            )
    from Vacuum_Polarization.main_script import (
            core_iteration,
            constant_lambda_iterations,
            walkback,
            full_script,
            )


