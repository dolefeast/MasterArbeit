from itertools import count
def main(
    m = 0,
    lambda_min = 16,
    lambda_max = 30,
    lambda_step = 0.105,
    n_iterations = None,
    verbose = 1,
    tol = 5e-5,
    max_nodes=5e8,
    plot: bool=True,
    save_solutions=False,
    read_solutions=True,
    coupling: float=1,
        ):

    from physics import Vacuum_Polarization

    if not save_solutions:
        print('Watch out! Solutions will NOT be saved')

    if plot:
        import matplotlib.pyplot as plt
        import matplotlib

        fig, (
                ax_rho,
                ax_A0_induced,
                ax_eigenvalue_array,
                )= plt.subplots(3)
        plot_rho = plot
        plot_A0_induced = plot
    else:
        ax_rho = None
        ax_A0_induced = None

    system = Vacuum_Polarization(
            lambda_value=lambda_min,
            m=m,
            read_solutions=read_solutions,
            coupling=coupling,
            )

    if plot:
        eigenvalue_array_array = []
        lambda_value_array = []
    if system.read_solutions:
        lambda_min = lambda_min + lambda_step
        system.lambda_value = lambda_min
    try:

        for lambda_value in count(lambda_min, lambda_step):
            print(f'Calculating solutions for lambda_value = {lambda_value}')
            system.update_eigenstates_script(
                    n_iterations=n_iterations,
                    verbose=verbose,
                    plot_rho=plot_rho,
                    ax_rho=ax_rho,
                    plot_A0_induced=plot_A0_induced,
                    ax_A0_induced=ax_A0_induced,
                    save_solutions=save_solutions,
                    tol=tol,
                    max_nodes=max_nodes,
                    )
            if system.broken:
                break
            system = Vacuum_Polarization(
                    lambda_value=lambda_value,
                    m=m,
                    read_solutions=False,
                    eigenvalue_array=system.eigenvalue_array,
                    eigenstate_array=system.eigenstate_array,
                    eigenstate_gradient_array=system.eigenstate_gradient_array,
                    A0_induced=system.A0_induced,
                    )
            system.lambda_value = lambda_value
            if plot:
                eigenvalue_array_array.append(system.eigenvalue_array)
                lambda_value_array.append(system.lambda_value)
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

    return system.lambda_value
