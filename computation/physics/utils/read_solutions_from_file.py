import numpy as np

from utils.float_to_str import  float_to_str

def read_solutions_from_file(self):
    m_string = float_to_str(self.m, self.sig_digs)
    a_string = float_to_str(self.a, self.sig_digs)
    E_string = float_to_str(self.E, self.sig_digs)

    directory = f"{self.read_solutions_dir}/{self.bcs}",
    # file_id = f"mass_{m_string}_lambda_{lambda_string}.txt"
    file_id = f"mass_{m_string}_a_{a_string}_E_{E_string}.txt"

    self.rho_array  = np.genfromtxt(
        f"{directory}/rho/{file_id}",
        dtype=float,
        delimiter="\n",
    )

    self.A0_induced = np.genfromtxt(
        f"{directory}/A0_induced/{file_id}",
        dtype=float,
        delimiter=",",
    )

    self.A0_field.value = -self.lambda_value * (self.z - 1/2) + self.A0_induced

    self.eigenstate_array =  list(
        np.genfromtxt(
        f"{directory}/eigenstate_array/{file_id}", delimiter=","
        )
    )
    self.eigenstate_gradient_array = list(
        np.genfromtxt(
        f"{directory}/eigenstate_gradient_array/{file_id}",
        delimiter=",",
        )
    )

    self.eigenvalue_array = list(
        np.genfromtxt(
        f"{directory}/eigenvalue_array/{file_id}", delimiter="\n"
        )
    )
