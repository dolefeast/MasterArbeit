class Vacuum_Polarization:

    from physics.utils.Klein_Gordon import (
            Klein_Gordon,
            dirichlet_boundary_conditions,
            neumann_boundary_conditions
            )

    from physics.init.init_method import (
            __init__,
            lambda_value
            )

    # To read possible existing solutions
    from physics.utils.read_solutions_from_file import read_solutions_from_file

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

    # To normalize the eigenstates
    from physics.calculate.calculate_rho import (
            eigenstate_rho,
            calculate_total_rho,
            )

