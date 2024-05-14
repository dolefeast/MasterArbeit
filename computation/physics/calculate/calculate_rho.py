from numpy import zeros_like, abs, pi

def eigenstate_rho(
        self,
        eigenvalue,
        eigenstate,
        ):
    return (
            (eigenvalue - self.e * self.A0_field.value) 
            * abs(eigenstate)**2
            )

def calculate_total_rho(
        self,
        ):
    """
    Calculates the corresponding charge density rho corresponding to self.eigenstate_array.
    Watch out! Assumes normalized eigenstates!
    """
    total_rho = zeros_like(self.z, dtype='float64')
    for eigenvalue, eigenstate in zip(
            self.eigenvalue_array,
            self.eigenstate_array
            ):
        # Calculate the [float] format of the charge density 
        # associated to the nth eigenstate
        rho_n = self.eigenstate_rho(eigenvalue, eigenstate)

        total_rho = total_rho + rho_n

    self.rho = 1/2 * total_rho + self.e ** 2 / pi * self.A0_field.value
