import numpy as np

from math_objects.perturbative_solutions import dirichlet_eigenstate, dirichlet_eigenstate_gradient
from math_objects.Klein_Gordon import *
from math_objects.fields import Vector_Potential


def __init__(self,
        lambda_value: float=None,
        A0_perturbation: [float]=None,
        m: float=1,
        e: float=1,
        n_points=500,
        scalar_name: str = "phi1",
        eigenvalue_array: [float]=None,
        eigenstate_array=None,
        eigenstate_gradient_array=None,
        N_mode_cutoff: int = 150,
        float_tol=1e-2,
        boundary_conditions: callable = dirichlet_boundary_conditions,
        read_solutions: bool=True,
        sig_digs: int=3,
        ):
    #Computation stuff:
    self.n_points = n_points
    self.float_tol = float_tol
    self.scalar_name = scalar_name
    self.scalar_name = scalar_name
    self.broken = 0 # When calculating eigenstates
                    # if an eigenvalue gets 
                    # repeated, or the solution
                    # does not converge, stop the
                    # calculation

    #Physics stuff:
    self.e = e
    self.z = np.linspace(0, 1, n_points)
    self.N_mode_cutoff = N_mode_cutoff
    self.sig_digs = sig_digs # Significative digits
    self.lambda_value = round(lambda_value, sig_digs)
    self.m = round(m, sig_digs)

    self.A0_field = Vector_Potential(
    n_points=self.n_points
            )
    self.boundary_conditions = self.dirichlet_boundary_conditions

    if eigenstate_array is None:
        try:
            # Check and read if corresponding file exists
            if read_solutions:
                self.read_solutions()
            else: 
                raise FileNotFoundError
        except FileNotFoundError:
            print('Tried to read the solutions for lambda = {self.lambda_value} and mass = {self.m} but could not find it.\n\tEigenstates will be generated assuming Dirichlet boundary conditions.')
            # File did not exist => generate guesses
            if not eigenvalue_array is None:
                self.eigenvalue_array = eigenvalue_array
            else: 
                n = np.arange(-N_mode_cutoff, N_mode_cutoff+1)
                # lambda=0 case is a super good enough approximation
                self.eigenvalue_array = np.sign(n) * np.sqrt( 
                        n**2 * np.pi**2
                        + self.m**2
                        )
                        
            # Guesses for the boundary value problem solution.
            print("Warning: No eigenstate array was given. It will be created assuming Dirichlet boundary conditions.")
            self.eigenstate_array = [
                    dirichlet_eigenstate(self.z, omega, self.m, self.lambda_value)
                    for omega in self.eigenvalue_array
                    ]
            self.eigenstate_gradient_array = [
                    dirichlet_eigenstate_gradient(self.z, omega, self.m, lambda_value)
                    for omega in self.eigenvalue_array
                    ]
    else:
        if eigenvalue_array is None:
            # if it is not given, retrieve it from the solutions array
            # this is outdated. i dont use bvp solutions anymore
            self.eigenvalue_array = [sol.p[0] for sol in eigenstate_array]
        else:  # If the eigenvalue_array is given
            # They must have same length.
            assert len(eigenvalue_array) == len(eigenstate_array)
            self.eigenvalue_array = eigenvalue_array
        # not anymore interested in the solve_bvp.solution format
        self.eigenstate_array = eigenstate_array 
        self.eigenstate_gradient_array = eigenstate_gradient_array 

    self.A0_base_value = - self.lambda_value/self.e * (self.z - 1/2)

    if not A0_perturbation is None:
        self.A0_perturbation = A0_perturbation
        self.A0_value = self.A0_base_value + self.A0_perturbation
    else: 
        self.A0_value = self.A0_base_value

    self.A0_field = Vector_Potential(
            value=self.A0_value, n_points=self.n_points
            )
