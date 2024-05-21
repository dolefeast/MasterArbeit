from itertools import count

def main(
        lambda_min=14.5,
        lambda_max=16,
        m=0,
        n_iterations=None,
        verbose=0,
        plot=False,
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
            sig_digs=sig_digs,
            coupling=coupling,
            )

    lambda_array_array = [
            np.linspace(lambda_min, lambda_max, lambda_step_count)
            for lambda_step_count in range(2, 15)
            ]

    for lambda_array in lambda_array_array:
        try:
            for lambda_value in lambda_array:
                if verbose:
                    print(f'Calculating convergence for lambda_value={system.lambda_value}')
                system.update_eigenstates_script(
                        n_iterations=n_iterations,
                        verbose=verbose,
                        plot_rho=False,
                        ax_rho=ax_rho,
                        plot_A0_induced=False,
                        ax_A0_induced=ax_A0_induced,
                        save_solutions=save_solutions,
                        tol=tol,
                        max_nodes=max_nodes,
                        )
                if system.broken:
                    system.plot_rho_A0(
                            True,
                            ax_rho,
                            True,
                            ax_A0_induced,
                            1,
                            )
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
        fig.savefig('/home/sanz/lambda_step_comparison.pdf')
        plt.show()

    try:
        raise exception
    except UnboundLocalError:
        pass

    return system.lambda_value
