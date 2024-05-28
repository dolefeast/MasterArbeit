from itertools import count

def main(
        m,
        a,
        lambda_min,
        lambda_step,
        n_iterations,
        verbose,
        plot,
        save_solutions,
        read_solutions,
        lambda_max=None,
        directory="",
        tol=1e-4,
        max_nodes=5e6,
        sig_digs=4,
        coupling=1,
        ):

    from physics import Vacuum_Polarization
    import numpy as np

    if plot:
        import matplotlib.pyplot as plt
        fig, (ax_rho, ax_A0_induced, ax_eigenvalue_array) = plt.subplots(3)
        eigenvalue_array_array = []
        lambda_value_array = []
    else:
        ax_rho = None
        ax_A0_induced = None
        ax_eigenvalue_array = None

    system = Vacuum_Polarization(
            lambda_value=lambda_min,
            a=a,
            m=m,
            read_solutions=read_solutions,
            read_solutions_dir=directory,
            save_solutions_dir=directory,
            )

    # Will be 0 if the file was not found
    if system.read_solutions:
        lambda_min = lambda_min + 2 * lambda_step
        system.lambda_value = lambda_min - lambda_step

    try:
        for lambda_value in count(lambda_min, lambda_step):
            if isinstance(lambda_max, (int, float)):
                if lambda_value>lambda_max:
                    break
            print('Calculating for lambda_value=', system.lambda_value)
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
            if plot:
                eigenvalue_array_array.append(system.eigenvalue_array)
                lambda_value_array.append(system.lambda_value)
            system.lambda_value = lambda_value
    # Whatever happens here, break the loop and plot whatever we have.
    except KeyboardInterrupt:
        pass
    except Exception as e:
        exception = e
        pass

    if plot:
        ax_eigenvalue_array.plot(
        lambda_value_array,
        eigenvalue_array_array,
        'b'
        )

        ax_rho.set_ylabel(r"$\rho(z)$")
        ax_A0_induced.set_ylabel(r"$A_0(z)$")
        ax_eigenvalue_array.set_xlabel(r"$\lambda$")
        ax_eigenvalue_array.set_ylabel(r"$\omega_n$")
        plt.show()

    try:
        raise exception
    except UnboundLocalError:
        pass

    return system.lambda_value - lambda_step
