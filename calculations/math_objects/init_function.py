import numpy as np

from math_objects.perturbative_solutions import (
    dirichlet_eigenstate,
    dirichlet_eigenstate_gradient,
    neumann_eigenstate,
    neumann_eigenstate_gradient,
)
from math_objects.Klein_Gordon import *
from math_objects.fields import Vector_Potential


def __init__(
    self,
    lambda_value: float,
    A0_induced: [float]=None,
    m: float=1,
    e: float=1,
    n_points=500,
    scalar_name: str="phi1",
    eigenvalue_array: [float]=None,
    eigenstate_array=None,
    eigenstate_gradient_array=None,
    N_mode_cutoff: int=150,
    float_tol=1e-2,
    bcs = 'dirichlet',
    read_solutions: bool=True,
    sig_digs: int=3,
):
    # Computation stuff:
    self.n_points = n_points
    self.float_tol = float_tol
    self.scalar_name = scalar_name
    self.scalar_name = scalar_name
    self.broken = 0  # When calculating eigenstates
    # if an eigenvalue gets
    # repeated, or the solution
    # does not converge, stop the
    # calculation

    # Physics stuff:
    self.e = e
    self.z = np.linspace(0, 1, n_points)
    self.N_mode_cutoff = N_mode_cutoff
    self.sig_digs = sig_digs  # Significative digits
    self.lambda_value = round(lambda_value, sig_digs)
    self.m = round(m, sig_digs)

    self.boundary_conditions = self.dirichlet_boundary_conditions
    self.bcs = bcs
    # Guesses for the boundary value problem solution.
    if bcs == 'dirichlet':
        self.perturbative_eigenstate = dirichlet_eigenstate
        self.perturbative_eigenstate_gradient = dirichlet_eigenstate_gradient
    elif bcs == 'neumann':
        self.perturbative_eigenstate = neumann_eigenstate
        self.perturbative_eigenstate_gradient = neumann_eigenstate_gradient
    else:
        self.bcs = 'dirichlet'
        print('No boundary conditions were given, assuming dirichlet bcs')

    if A0_induced is None:
        self.A0_induced = np.zeros(self.n_points)
    else:
        self.A0_induced = A0_induced

    self.A0_base_value = -self.lambda_value * (self.z - 1/2)
    self.A0_field = Vector_Potential(
    n_points=self.n_points,
    value = -self.A0_base_value + self.A0_induced
            )

    self.eigenstate_array = eigenstate_array
    self.eigenstate_gradient_array = eigenstate_gradient_array
    self.eigenvalue_array = eigenvalue_array
    self.read_solutions = read_solutions

    self.generate_eigenstates()

