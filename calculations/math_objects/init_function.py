import numpy as np
from math_objects.Klein_Gordon import *
from math_objects.fields import Vector_Potential

def __init__(self,
        lambda_value: float=None,
        A0_modification: [float]=None,
        m: float=1,
        e: float=1,
        n_points=500,
        scalar_name: str = "phi1",
        eigenvalue_guess_array: [float]=None,
        eigenstate_array=None,
        N_mode_cutoff: int = 150,
        float_tol=1e-2,
        boundary_conditions: callable = dirichlet_boundary_conditions
        ):
    #Computation stuff:
    self.n_points = n_points
    self.float_tol = float_tol

    #Physics stuff:
    self.m = m
    self.e = e
    self.z = np.linspace(0, 1, n_points)
    self.scalar_name = scalar_name
    if not eigenstate_array is None: 
        if eigenvalue_guess_array is None:
            # if it is not given, calculate it from the solutions array
            self.eigenvalue_guess_array = [sol.p[0] for sol in eigenstate_array]
        else:  # If the eigenvalue_guess_array is given
            # They must have same length.
            assert len(eigenvalue_guess_array) == len(eigenstate_array)
        # not anymore interested in the solve_bvp.solution format
        self.eigenstate_array = [sol.y[0] for sol in eigenstate_array] 
        self.eigenstate_gradient_array = [sol.y[1] for sol in eigenstate_array] 
    else:
        if not eigenvalue_guess_array is None:
            self.eigenvalue_guess_array = eigenvalue_guess_array
        else: 
            n = np.arange(-N_mode_cutoff, N_mode_cutoff+1)
            # lambda=0 case is a good enough approximation
            self.eigenvalue_guess_array = np.sign(n) * np.sqrt( 
                    n**2 * np.pi**2
                    + self.m**2
                    )
        # Guesses for the boundary value problem solution.
        self.eigenstate_array = [1/omega * np.sin(omega*self.z) 
                for omega in self.eigenvalue_guess_array]
        print("Warning: Since no eigenstate array was given, creating it assuming Dirichlet boundary conditions")
        self.eigenstate_gradient_array = [np.cos(omega*self.z) 
                for omega in self.eigenvalue_guess_array]

    self.boundary_conditions = self.dirichlet_boundary_conditions

    self.z = np.linspace(0, 1, n_points)

    self.lambda_value = lambda_value
    self.A0_base_value = - self.lambda_value/self.e * (self.z - 1/2)

    if not A0_modification is None:
        if np.abs(A0_modification[0]) < float_tol or  np.abs(A0_modification[-1]) < float_tol: 
            raise AssertionError("The electric potential is not 0 in the boundaries of the problem.")

        # This can be done better, I don't think there's a need to, right now.
        self.A0_modification = A0_modification

        self.A0_value = self.A0_base_value + self.A0_modification
    else: 
        self.A0_value = self.A0_base_value

    self.A0_field = Vector_Potential(value=self.A0_value, n_points=self.n_points
            )
