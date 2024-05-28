from numpy import savetxt
from utils.float_to_str import float_to_str

def save_solutions(self):
    lambda_string = float_to_str(self.lambda_value, self.sig_digs)
    a_string = float_to_str(self.a, self.sig_digs)
    m_string = float_to_str(self.m, self.sig_digs)

    file_id = f"mass_{m_string}_a_{a_string}_lambda_{lambda_string}.txt"

    if self.save_solutions_dir != "" and self.save_solutions_dir[0]!="/": 
        save_solutions_dir = "/" + self.save_solutions_dir 

    main_directory = f"saved_solutions{save_solutions_dir}/{self.bcs}"

    print(f"Saving results under {directory + '/.../' + file_id}")
    # Saving the eigenvalues
    to_csv = self.eigenvalue_array
    savetxt(
            f"{main_directory}/eigenvalue_array/{file_id}",
            to_csv, 
            delimiter=","
            )

    # Saving the eigenstates
    to_csv = self.eigenstate_array
    savetxt(
            f"{main_directory}/eigenstate_array/{file_id}",
        to_csv,
        delimiter=",",
    )

    # Saving the eigenstates gradient
    to_csv = self.eigenstate_gradient_array
    savetxt(
            f"{main_directory}/eigenstate_gradient_array/{file_id}",
        to_csv,
        delimiter=",",
    )

    # Saving the perturbation in the field
    to_csv = self.A0_induced
    # to_csv = self.A0_field.value + self.lambda_value * (self.z - 1 / 2)
    savetxt(
            f"{main_directory}/A0_induced/{file_id}",
             to_csv, delimiter=",")

    # Saving the resulting charge density
    to_csv = self.rho
    savetxt(
            f"{main_directory}/rho/{file_id}",
            to_csv,
            delimiter=","
            )
