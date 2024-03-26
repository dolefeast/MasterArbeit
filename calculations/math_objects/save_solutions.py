import numpy as np
from scripts.float_to_str import float_to_str
def save_solutions(self):
        lambda_string = float_to_str(self.lambda_value, self.sig_digs)
        m_string = float_to_str(self.m, self.sig_digs)
        file_id = f'lambda_{lambda_string}_mass_{m_string}.txt'
        print(f"Saving results under {file_id}")

        # Saving the eigenvalues 
        to_csv = self.eigenvalue_array
        np.savetxt(f'saved_solutions/dirichlet/eigenvalue/{file_id}', to_csv, delimiter= ",")

        # Saving the eigenstates 
        to_csv = self.eigenstate_array
        np.savetxt(f'saved_solutions/dirichlet/normalized_eigenstate/{file_id}', to_csv, delimiter= ",")

        # Saving the eigenstates gradient
        to_csv = self.eigenstate_gradient_array
        np.savetxt(f'saved_solutions/dirichlet/normalized_eigenstate_gradient/{file_id}', to_csv, delimiter= ",")

        # Saving the perturbation in the field
        to_csv = self.A0_field.value + self.lambda_value * (self.z - 1/2)
        np.savetxt(f'saved_solutions/dirichlet/A0_field/{file_id}', to_csv, delimiter= ",")

        # Saving the resulting charge density
        to_csv = self.rho_array
        np.savetxt(f'saved_solutions/dirichlet/rho/{file_id}', to_csv, delimiter= ",")
