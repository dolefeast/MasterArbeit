from numpy import savetxt
from utils.float_to_str import float_to_str

def save_solutions(self):
    lambda_string = float_to_str(self.lambda_value, self.sig_digs)
    m_string = float_to_str(self.m, self.sig_digs)

    file_id = f"mass_{m_string}_lambda_{lambda_string}.txt"
    print(f"Saving results under {file_id}")

    # Saving the eigenvalues
    to_csv = self.eigenvalue_array
    savetxt(f"saved_solutions/{self.bcs}/eigenvalue_array/{file_id}", to_csv, delimiter=",")

    # Saving the eigenstates
    to_csv = self.eigenstate_array
    savetxt(
        f"saved_solutions/{self.bcs}/eigenstate_array/{file_id}",
        to_csv,
        delimiter=",",
    )

    # Saving the eigenstates gradient
    to_csv = self.eigenstate_gradient_array
    savetxt(
        f"saved_solutions/{self.bcs}/eigenstate_gradient_array/{file_id}",
        to_csv,
        delimiter=",",
    )

    # Saving the perturbation in the field
    to_csv = self.A0_field.value + self.lambda_value * (self.z - 1 / 2)
    savetxt(f"saved_solutions/{self.bcs}/A0_induced/{file_id}", to_csv, delimiter=",")

    # Saving the resulting charge density
    to_csv = self.rho
    savetxt(f"saved_solutions/{self.bcs}/rho/{file_id}", to_csv, delimiter=",")
