from itertools import count
def main(
    m = 0,
    lambda_min = 16,
    lambda_step = 0.01,
    n_iterations = None,
    verbose = 1,
    tol = 5e-5,
    max_nodes=5e8,
    plot_rho=True,
    plot_A0_induced=True,
    save_solutions=False,
    read_solutions=True
        ):

    from physics import Vacuum_Polarization

    if not save_solutions:
        print('Watch out! Solutions will NOT be saved')

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
        system.lambda_vlue = lambda_min
    try:

        for lambda_value in count(lambda_min, lambda_step):
            print(f'Calculating solutions for lambda_value = {system.lambda_value}')
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

    return system.lambda_value
