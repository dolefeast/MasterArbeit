from itertools import count
def main(
        m=0,
        lambda_min=16,
        lambda_max=None,
        lambda_step=0.1,
        n_iterations=None,
        verbose=1,
        plot=True,
        save_solutions=False,
        read_solutions=True,
        tol=1e-4,
        max_nodes=5e6,
        sig_digs=4,
        coupling=1,
        ):

    from physics import Vacuum_Polarization
    import numpy as np

    if plot:
        import matplotlib.pyplot as plt
        fig, (ax_rho, ax_A0_induced) = plt.subplots(2)
    else:
        ax_rho = None
        ax_A0_induced = None

    system = Vacuum_Polarization(
            lambda_value=lambda_min,
            m=m,
            read_solutions=read_solutions,
            sig_digs=4,
            )


    if system.read_solutions:
        system.lambda_value = lambda_min + lambda_step
        lambda_min = lambda_min + 2*lambda_step

    if isinstance(lambda_max, float):
        lambda_array = np.arange(lambda_min, lambda_max, lambda_step)
    else:
        lambda_array = count(lambda_min, lambda_step)

    try:
        for lambda_value in lambda_array:
            if verbose:
                print(f'Calculating convergence for lambda_value={system.lambda_value}')
            system.update_eigenstates_script(
                    n_iterations=n_iterations,
                    verbose=verbose,
                    plot_rho=plot,
                    ax_rho=ax_rho,
                    plot_A0_induced=plot,
                    ax_A0_induced=ax_A0_induced,
                    save_solutions=save_solutions,
                    tol=tol,
                    max_nodes=max_nodes,
                    )
            if system.broken:
                break
            system.lambda_value = lambda_value
    # Whatever happens here, break the loop and plot whatever we have.
    except KeyboardInterrupt:
        pass
    except Exception as e:
        exception = e
        pass

    if plot:
        ax_rho.set_ylabel(r"$\rho(z)$")
        ax_A0_induced.set_ylabel(r"$A_0(z)$")
        plt.show()

    try:
        raise exception
    except UnboundLocalError:
        pass

    return system.lambda_value
