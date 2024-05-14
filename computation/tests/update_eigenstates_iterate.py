def main(
        m,
        lambda_value,
        n_iterations,
        verbose,
        plot_rho,
        plot_A0_induced,
        ):

    from physics import Vacuum_Polarization

    if plot_rho or plot_A0_induced:
        import matplotlib.pyplot as plt
        fig, (ax_rho, ax_A0_induced) = plt.subplots(2)

    system = Vacuum_Polarization(
            lambda_value=lambda_value,
            m=m,
            )

    system.update_eigenstates_iterate(
            n_iterations,
            verbose,
            plot_rho,
            ax_rho,
            plot_A0_induced,
            ax_A0_induced,
            )

    if plot_rho:
        ax_rho.set_ylabel(r"$\rho(z)$")

    if plot_A0_induced:
        ax_A0_induced.set_ylabel(r"$A_0(z)$")

    if plot_rho or plot_A0_induced:
        plt.show()
