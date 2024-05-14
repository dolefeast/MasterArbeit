def main(
        m,
        lambda_value,
        n_iterations,
        verbose,
        plot_rho,
        plot_A0_induced,
        save_solutions,
        tol=1e-3,
        ):

    from physics import Vacuum_Polarization

    if plot_rho or plot_A0_induced:
        import matplotlib.pyplot as plt
        fig, (ax_rho, ax_A0_induced) = plt.subplots(2)
    else:
        ax_rho = None
        ax_A0_induced = None

    system = Vacuum_Polarization(
            lambda_value=lambda_value,
            m=m,
            )

    system.update_eigenstates_script(
            n_iterations=n_iterations,
            verbose=verbose,
            plot_rho=plot_rho,
            ax_rho=ax_rho,
            plot_A0_induced=plot_A0_induced,
            ax_A0_induced=ax_A0_induced,
            save_solutions=save_solutions,
            tol=tol,
            )

    if plot_rho:
        ax_rho.set_ylabel(r"$\rho(z)$")

    if plot_A0_induced:
        ax_A0_induced.set_ylabel(r"$A_0(z)$")

    if ( plot_rho or plot_A0_induced ) and not system.broken:
        plt.show()
    elif system.broken:
        print('The calculation broke at some point. No output plot is generated')
