import numpy as np

from utils.float_to_str import  float_to_str

def read_solutions_from_file(self):
    m_string = float_to_str(self.m, self.sig_digs)
    lambda_string = float_to_str(self.lambda_value, self.sig_digs)
    a_string = float_to_str(self.a, self.sig_digs)
    file_id = f"mass_{m_string}_a_{a_string}_lambda_{lambda_string}.txt"

    if self.read_solutions_dir != "" and self.read_solutions_dir[0]!="/": 
        read_solutions_dir = "/" + self.read_solutions_dir 
    else:
        read_solutions_dir = self.read_solutions_dir

    self.rho_array  = np.genfromtxt(
        f"saved_solutions{read_solutions_dir}/{self.bcs}/rho/{file_id}",
        dtype=float,
        delimiter="\n",
    )

    self.A0_induced = np.genfromtxt(
        f"saved_solutions{read_solutions_dir}/{self.bcs}/A0_induced/{file_id}",
        dtype=float,
        delimiter=",",
    )

    self.A0_field.value = -self.lambda_value * (self.z - 1/2) + self.A0_induced

    self.eigenstate_array =  list(
        np.genfromtxt(
        f"saved_solutions{read_solutions_dir}/{self.bcs}/eigenstate_array/{file_id}", delimiter=","
        )
    )
    self.eigenstate_gradient_array = list(
        np.genfromtxt(
        f"saved_solutions{read_solutions_dir}/{self.bcs}/eigenstate_gradient_array/{file_id}",
        delimiter=",",
        )
    )

    self.eigenvalue_array = list(
        np.genfromtxt(
        f"saved_solutions{read_solutions_dir}/{self.bcs}/eigenvalue_array/{file_id}", delimiter="\n"
        )
    )
