import numpy as np

from scripts.float_to_str import  float_to_str
from scripts.antisymmetry import  antisymmetry, antisymmetry_eigenstates
from math_objects.fields import Vector_Potential


def read_solutions(self):
    m_string = float_to_str(self.m, self.sig_digs)
    lambda_string = float_to_str(self.lambda_value, self.sig_digs)
    file_id = f"lambda_{lambda_string}_mass_{m_string}.txt"

    self.rho_array  = antisymmetry(
            np.genfromtxt(
            f"saved_solutions/dirichlet/rho/{file_id}",
            dtype=float, 
            delimiter="\n",
            )
        )    
    A0_modification = np.genfromtxt(
        f"saved_solutions/dirichlet/A0_field/{file_id}",
        dtype=float,
        delimiter="\n",
    )



    self.A0_field.value = -self.lambda_value * (self.z - 1/2) + A0_modification


    self.eigenstate_array = np.genfromtxt(
        f"saved_solutions/dirichlet/normalized_eigenstate/{file_id}", delimiter=","
    )

    self.eigenstate_gradient_array = np.genfromtxt(
        f"saved_solutions/dirichlet/normalized_eigenstate_gradient/{file_id}",
        delimiter=",",
    )
    self.eigenvalue_array = np.genfromtxt(
        f"saved_solutions/dirichlet/eigenvalue/{file_id}", delimiter="\n"
    )
