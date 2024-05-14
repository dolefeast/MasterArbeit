class Vacuum_Polarization:

    from physics.utils.Klein_Gordon import (
            Klein_Gordon,
            dirichlet_boundary_conditions,
            neumann_boundary_conditions
            )

    from physics.init.init_method import (
            __init__,
            lambda_value,
            E_induced,
            )

    # To read possible existing solutions
    from physics.utils.read_solutions_from_file import read_solutions_from_file
    from physics.utils.save_solutions import save_solutions
    from physics.utils.plot_rho_A0 import plot_rho_A0

    from physics.filtering.convolve import (
            convolve,
            convolve_twice,
            )

    from physics.filtering.total_filtering import (
            total_filtering_dirichlet
            )

    from physics.filtering.filter_scripts import (
            extend_signal,
            remove_neighbourhood,
            remove_and_interpolate,
            return_to_0_1,
            )

    from physics.init.define_eigenstates import (
            # Checks whether the eigenstates are given, should be read or calculated
            define_eigenstates,
            # Calculates perturbative guesses for the eigenstates
            calculate_perturbative_eigenstates
            )

    # Calculates the eigenstates of the Klein-Gordon equation
    from physics.calculate.calculate_eigenstates import calculate_eigenstates

    # To normalize the eigenstates
    from physics.calculate.normalize_eigenstates import (
            eigenstate_norm,
            normalize_eigenstates
            )

    # To calculate the total charge density rho
    from physics.calculate.calculate_rho import (
            eigenstate_rho,
            calculate_total_rho,
            )

    from physics.calculate.filter_rho import (
            filter_rho,
            )

    from physics.calculate.update_eigenstates import (
            update_eigenstates,
            )

    from physics.calculate.update_eigenstates_iterate import (
            update_eigenstates_iterate,
            )

    from physics.calculate.update_eigenstates_script import (
            update_eigenstates_script,
            )

    from physics.calculate.update_eigenstates_converge import (
            update_eigenstates_converge,
            )

    from physics.calculate.calculate_A0_induced import (
            calculate_A0_induced,
            )

