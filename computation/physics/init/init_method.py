import numpy as np
from itertools import cycle

from physics.utils.Fields import Field
from physics.utils.perturbative_solutions import (
    dirichlet_eigenstate,
    dirichlet_eigenstate_gradient,
    neumann_eigenstate,
    neumann_eigenstate_gradient,
)

def __init__(
        self,
        lambda_value: float,
        a: float=1,
        m: float=0,
        e: float=1,
        eigenvalue_array: [float]=None,
        eigenstate_array: [float]=None,
        eigenstate_gradient_array: [float]=None,
        A0_induced: [float]=0,
        n_points: int=400,
        N_mode_cutoff: int=49,
        bcs: str='dirichlet',
        read_solutions: bool=False,
        sig_digs: int=3,
        float_tol: int=1e-2,
        read_solutions_dir:str="",
        save_solutions_dir:str="",
        ambjorn: bool=False,
        ):
    # Computation 
    self.n_points = n_points
    self.N_mode_cutoff = N_mode_cutoff
    self.broken = 0 # Tracks whether the calculation broke down at any point
    self.float_tol = float_tol # To track if the calculation broke down at any point

    self.e = e
    self.a = a

    self._lambda_value = round(lambda_value, sig_digs)
    self.m = round(m, sig_digs)
    self.z = np.linspace(0, 1, n_points)
    self.sig_digs = sig_digs
    # If True, calculates the charge density as Ambjorn and Wolfram do
    self.ambjorn = ambjorn

    self.read_solutions_dir = read_solutions_dir 
    self.save_solutions_dir = save_solutions_dir

    self.bcs = bcs
    if bcs == 'dirichlet':
        self.perturbative_eigenstate = dirichlet_eigenstate
        self.perturbative_eigenstate_gradient = dirichlet_eigenstate_gradient
        self.boundary_conditions = self.dirichlet_boundary_conditions
    elif bcs == 'neumann':
        self.perturbative_eigenstate = neumann_eigenstate
        self.perturbative_eigenstate_gradient = neumann_eigenstate_gradient
        self.boundary_conditions = self.neumann_boundary_conditions

    if isinstance(A0_induced, (int, float)):
        self.A0_induced = np.array(
                [A0_induced]*self.n_points
                )
    else:
        self.A0_induced = np.array(A0_induced)

    self._E_induced = None

    # Need to define Field because it needs to be callable
    # for the klein gordon equation
    self.A0_field = Field(
            n_points=self.n_points,
            value = (
                - self.lambda_value 
                * (self.z - 1/2)
                + self.a * self.A0_induced
                )
            )

    # These may or may not be None
    self.eigenstate_array = eigenstate_array
    self.eigenstate_gradient_array = eigenstate_gradient_array
    self.eigenvalue_array = eigenvalue_array
    self.read_solutions = read_solutions

    # Define the eigenstate arrays
    self.define_eigenstates()

@property
def lambda_value(self):
    return self._lambda_value

@lambda_value.setter
def lambda_value(
        self, 
        new_lambda_value
        ):
    self._lambda_value = new_lambda_value
    self.A0_field = Field(
            n_points=self.n_points,
            value = (
            - self.lambda_value 
            * (self.z - 1/2)
            + self.a * self.A0_induced
        )
    )

@property
def E_induced(
        self,
        ):
    if self._E_induced is None:
        return -np.diff(self.A0_induced)
    return self._E_induced
    

