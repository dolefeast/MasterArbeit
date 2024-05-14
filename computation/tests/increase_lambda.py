from itertools import count
def main(
        m,
        lambda_min,
        lambda_step,
        n_iterations,
        verbose,
        plot_rho,
        plot_A0_induced,
        save_solutions,
        read_solutions,
        tol=1e-4,
        ):

    from physics import Vacuum_Polarization

    if plot_rho or plot_A0_induced:
        import matplotlib.pyplot as plt
        fig, (ax_rho, ax_A0_induced) = plt.subplots(2)
    else:
        ax_rho = None
        ax_A0_induced = None

    system = Vacuum_Polarization(
            lambda_value=lambda_min,
            m=m,
            read_solutions=read_solutions
            )

    if system.read_solutions:
        lambda_min = lambda_min + lambda_step
        system.lambda_value = lambda_min
    try:
        system.A0_induced *= 0
        for lambda_value in count(lambda_min, lambda_step):
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
            if system.broken:
                break
            print('Increasing lambda_value to', lambda_value)
            system.lambda_value = lambda_value
    # Whatever happens here, break the loop and plot whatever we have.
    except KeyboardInterrupt:
        pass
    except Exception as e:
        exception = e
        pass

    if plot_rho:
        ax_rho.set_ylabel(r"$\rho(z)$")

    if plot_A0_induced:
        ax_A0_induced.set_ylabel(r"$A_0(z)$")

    if  plot_rho or plot_A0_induced :
        plt.show()

    try:
        raise exception
    except UnboundLocalError:
        pass
