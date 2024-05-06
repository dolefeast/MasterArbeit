import numpy as np

from utils.float_to_str import  float_to_str

def read_solutions_from_file(self):
    m_string = float_to_str(self.m, self.sig_digs)
    lambda_string = float_to_str(self.lambda_value, self.sig_digs)
    file_id = f"lambda_{lambda_string}_mass_{m_string}.txt"

    self.rho_array  = np.genfromtxt(
        f"saved_solutions/dirichlet/rho_array/{file_id}",
        dtype=float,
        delimiter="\n",
    )

    self.A0_induced = np.genfromtxt(
        f"saved_solutions/dirichlet/A0_induced/{file_id}",
        dtype=float,
        delimiter=",",
    )

    self.A0_field.value = -self.lambda_value * (self.z - 1/2) + self.A0_induced

    self.eigenstate_array =  list(
        np.genfromtxt(
        f"saved_solutions/dirichlet/eigenstate_array/{file_id}", delimiter=","
        )
    )
    self.eigenstate_gradient_array = list(
        np.genfromtxt(
        f"saved_solutions/dirichlet/eigenstate_gradient_array/{file_id}",
        delimiter=",",
    )
    )

    self.eigenvalue_array = list(
        np.genfromtxt(
        f"saved_solutions/dirichlet/eigenvalue_array/{file_id}", delimiter="\n"
        )
    )
