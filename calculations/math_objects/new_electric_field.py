import scipy as sp

def new_electric_field(self, eigenstate_array):
    """
    Calculates the corresponding electric field to a certain charge density
    Parameters:
        z: float. The z value at which the electric field is to be calculated
        eigenstate_array: SolutionArray. A given solution array to calculate_charge_density
    Returns:
        electric_field: float, the electric field at z
    """

    # Calculate this solving an ODE, not integration
    print("I am again here, hello!")
    charge_density_at_z = self.calculate_total_charge_density(eigenstate_array)

    modified_electric = sp.integrate.solve_ivp(
        lambda z, y: charge_density_at_z(z),
        y0=[0],
        t_span=(self.z[0], self.z[-1]),
        t_eval=self.z,
    )
    print("The ODE was completed!")

    return modified_electric

