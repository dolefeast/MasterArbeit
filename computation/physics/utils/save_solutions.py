from numpy import savetxt
from utils.float_to_str import float_to_str

def save_solutions(self,):
    """
    Saves the solutions to a certain string of the form 
    mass_0_0_a_1_0_E_5_0.txt 
    """

    E_string = float_to_str(self.E, self.sig_digs)
    a_string = float_to_str(self.a, self.sig_digs)
    m_string = float_to_str(self.m, self.sig_digs)

    directory = f"{self.save_solutions_dir}/{self.bcs}"
    file_id = f"mass_{m_string}_a_{a_string}_E_{E_string}.txt"


    print(f"Saving results under {directory + '/.../' + file_id}")
    # Saving the eigenvalues
    to_csv = self.eigenvalue_array
    savetxt(
            f"{directory}/eigenvalue_array/{file_id}", 
            to_csv,
            delimiter=","
            )

    # Saving the eigenstates
    to_csv = self.eigenstate_array
    savetxt(
        f"{directory}/eigenstate_array/{file_id}",
        to_csv,
        delimiter=",",
    )

    # Saving the eigenstates gradient
    to_csv = self.eigenstate_gradient_array
    savetxt(
        f"{directory}/eigenstate_gradient_array/{file_id}",
        to_csv,
        delimiter=",",
    )

    # Saving the perturbation in the field
    to_csv = self.A0_field.value + self.lambda_value * (self.z - 1 / 2)
    savetxt(f"{directory}/A0_induced/{file_id}", to_csv, delimiter=",")

    # Saving the resulting charge density
    to_csv = self.rho
    savetxt(f"{directory}/rho/{file_id}", to_csv, delimiter=",")
