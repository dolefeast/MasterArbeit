def main(
        m,
        lambda_value,
        ):
    from physics import Vacuum_Polarization
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(16, 9))

    system = Vacuum_Polarization(
            lambda_value=lambda_value,
            m=m,
            )
    mode = 50

    # KG solutions without normalization
    system.calculate_eigenstates()

    # Normalizing states
    system.normalize_eigenstates()

    # Calculate total charge density
    rho = system.calculate_total_rho()
    ax.plot(
            rho
            )

    rho_convolved = system.convolve(
            system.z,
            rho
            )

    ax.plot(
            rho_convolved
            )


    plt.show()
