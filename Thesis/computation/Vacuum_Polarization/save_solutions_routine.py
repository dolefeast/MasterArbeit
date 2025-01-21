import numpy as np
import os

def save_solutions(
        self,
        sig_digs=3,
        ):
    """
    Saves the calculated quantities to {directory}/{boundary_conditions}/{quantity}/mass_{mass}_a_{a}_lambda_value_{lambda_value}.csv
    Parameters:
        solution_family: A dictionary with keys() = ["eigenvalues", "eigenstates", "eigenstate_gradients"]
        directory: In case a further directory should be considered, e.g. if Ambjorn technique is used
    Returns None
    """

    solution_family = {
            "eigenvalues":self.eigenvalues,
            "eigenstates":self.eigenstates,
            "A0_induced":self.A0_induced(self.z),
            "rho":self.rho,
            }

    lambda_string = self.float_to_str(self.lambda_value, sig_digs=sig_digs)
    a_string = self.float_to_str(self.a, sig_digs=sig_digs)
    m_string = self.float_to_str(self.m, sig_digs=sig_digs)

    file_id = f"mass_{m_string}_a_{a_string}_lambda_{lambda_string}.txt"

    if self.directory != "":
        directory = "/" + self.directory

    root_directory = f"saved_solutions{directory}/{self.bcs}"
    print(f"Saving results under {root_directory}/.../{file_id}...")
    
    try:
        for key, value in solution_family.items():
            np.savetxt(
                    f"{root_directory}/{key}/{file_id}",
                value,
                delimiter=",",
            )
    except FileNotFoundError as e:
        print(e)
        create = input(f"\nCreate directory {directory[1:]}?[y/n]... ")
        if create == "y":
            for key in solution_family.keys():
                os.makedirs(root_directory + "/" + key)
            self.save_solutions()
        elif create == "n": 
            rename = input(f"If {directory} was a typo, enter the correct name...")
            if rename != "":
                self.save_solutions()
