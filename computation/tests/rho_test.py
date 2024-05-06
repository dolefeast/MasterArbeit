def main(
        lambda_min,
        ):
    from physics import Vacuum_Polarization
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    system = Vacuum_Polarization(lambda_min)
    mode = 50

    # KG solutions without normalization
    system.calculate_eigenstates()

    # Normalizing states
    system.normalize_eigenstates()

    eigenstate = system.eigenstate_array[mode]
    eigenvalue = system.eigenvalue_array[mode]

    ax.plot(
            system.eigenstate_rho(eigenvalue, eigenstate)
            )
    # Calculate total charge density
    ax.plot(
        system.calculate_total_rho()
            )

    plt.show()
