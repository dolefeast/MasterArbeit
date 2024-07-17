from math import pi
def update_eigenstates(
        self,
        smoothing: bool=True,
        max_nodes=5e6,
        ):
    """
    Given a state (A_0, phi_n)^kappa of the system, calculate the state kappa+1. 
    """
    # First calculate the (not normalized) KG solutions

    self.calculate_eigenstates(
            max_nodes=max_nodes,
            )

    # Normalizing the modes
    self.normalize_eigenstates()

    # Calculate total charge density
    self.calculate_total_rho()
    
    # Filter the charge density
    if smoothing: # ??
        self.filter_rho()

    if self.ambjorn:
        self.rho -=  self.e**2 / pi * self.A0_field.value
    
    # Calculate the corresponding A0
    self.calculate_A0_induced()
