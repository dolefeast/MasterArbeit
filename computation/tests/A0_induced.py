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
    alpha = 0.5

    # KG solutions without normalization
    system.calculate_eigenstates()

    # Normalizing states
    system.normalize_eigenstates()

    # Calculate total charge density
    rho = system.calculate_total_rho()

    z, rho = system.total_filtering_dirichlet(system.z, system.rho)

    ax.plot(system.z, rho, 
            alpha=alpha)

    system.calculate_A0_induced()
    ax.plot(system.z, system.A0_induced)

    plt.show()
