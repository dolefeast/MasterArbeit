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
    if smoothing:
        self.normalize_eigenstates()

    # Calculate total charge density
    self.calculate_total_rho()
    
    # Filter the charge density
    self.filter_rho()
    
    # Calculate the corresponding A0
    self.calculate_A0_induced()
