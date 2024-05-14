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
    alpha_unfiltered = 0.5
    rho = system.calculate_total_rho()

    # First extend
    z, filtered_rho = system.extend_signal(system.z, rho)

    ax.plot(
            z,
            filtered_rho,
            '--',
            alpha=alpha_unfiltered
            )

    # Then filter
    z, filtered_rho = system.convolve_twice(
            z,
            filtered_rho
            )

    ax.plot(
            z,
            filtered_rho,
            '--',
            alpha=alpha_unfiltered
            )

    # Take out 0 1 noisy neighbourhoods
    z, filtered_rho = system.remove_and_interpolate(
            z,
            filtered_rho
            )

    ax.plot(
            z,
            filtered_rho,
            '--',
            )
    # Return to z \in [0, 1]
    z, filtered_rho = system.return_to_0_1(
            z,
            filtered_rho
            )

    ax.plot(
            z,
            filtered_rho,
            label='step by step'
            )

    z, rho = system.total_filtering_dirichlet(
            z,
            rho
            )

    ax.plot(
            z,
            rho,
            '-.',
            label='doing the whole thing at once'
            )

    ax.legend(loc='best')
    plt.show()
