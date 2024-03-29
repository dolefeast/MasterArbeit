class Vacuum_Solution:
    """
    Parameters:
        lambda_value: float = None    #The the coupling constant of the electric potential
        A0_modification: np.array([float])xn_points = None  #The (external) electric potential in the case that is not a linear field
        m: float = 1      #The mass of the scalar field
        e: float = 1    #The charge of the scalar field
        n_points: int = 100         #The nodes to be considered
        scalar_name: str = 'phi'    #The name of the field, for plotting reasons
                                    #is subject to.
        eigenvalue_array: [float] # Array of eigenvalues for the bvp guess
        eigenstate_array: [solve_bvp.solutions] # Array of solutions from a previous calculation
        N_mode_cutoff: int=150 # If both eigenstate_array and eigenvalue_array are None, calculate 2 N_mode_cutoff modes (N_mode_cutoff positive and N_mode_cutoff negative energy solution)
        float_tol: float=1e-3 # The tolerance with which two floats are considered the same
    Properties:
        eigenstate_array: The eigenstate_array to a certain electric potential
    -------------
    Calculates the solution to the coupled differential equations
        (d²/dz² + (omega_n - eA_0(z))² - m²) phi = 0
        d²/dz² A = - rho
    In the context of a quantum charged scalar field phi and a classical electromagnetic potential (A_0, A_1) = (A_0, 0) under a certain gauge choice.

    The system to be considered. Takes as input the lamda value of the
    electric potential
        A_0(z) = - lambda_value * (z-1/2)
    and an array of eigenvalue omega_n guesses. These should be solutions to
    easier cases of the problem we are considering. For lamda small:
        omega_n = sign(n) * np.sqrt(n² pi² + m²)


    """

    from math_objects.init_function import __init__

    from math_objects.Klein_Gordon import (
        Klein_Gordon,
        dirichlet_boundary_conditions,
        neumann_boundary_conditions,
    )
    from math_objects.read_solutions import read_solutions
    from math_objects.save_solutions import save_solutions
    from math_objects.calculate_eigenstates import calculate_eigenstates
    from math_objects.charge_density import normalized_eigenstate, calculate_total_rho
    from math_objects.modify_A0 import modify_A0
    from math_objects.update_eigenstates import (
        update_eigenstates,
        update_eigenstates_iteration,
        update_eigenstates_until_convergence,
        update_eigenstates_script,
    )
