def main(
    m = 0,
    lambda_min = 16,
    inverse_lambda_step_min = 80,
    inverse_lambda_step_max = 120,
    inverse_lambda_step_count = 10,
    n_iterations = None,
    verbose = 1,
    tol = 5e-5,
    max_nodes=5e8,
    plot_rho=True,
    plot_A0_induced=True,
    save_solutions=False,
    read_solutions=True,
        ):

    import matplotlib.pyplot as plt
    from numpy import linspace

    from physics import Vacuum_Polarization
    from tests.increase_lambda import main as increase_lambda_main

    fig, ax = plt.subplots()

    inverse_lambda_step_array = linspace(
        inverse_lambda_step_min,
        inverse_lambda_step_max,
        inverse_lambda_step_count,
            )

    lambda_critical_array = []

    try:
        for inverse_lambda_step in inverse_lambda_step_array:
            lambda_critical = increase_lambda_main(
                m=m,
                lambda_min=lambda_min,
                lambda_step=inverse_lambda_step**-1,
                n_iterations=n_iterations,
                verbose=verbose,
                plot_rho=False,
                plot_A0_induced=False,
                save_solutions=False,
                read_solutions=True,
                tol=1e-4,
                max_nodes=5e6,
                    )
            lambda_critical_array.append(lambda_critical)
    except Exception as e:
        exception = e
        print(exception)
        pass
    except KeyboardInterrupt:
        pass

    ax.plot(
            inverse_lambda_step_array,
            lambda_critical_array,
            'x',
                        )

    ax.set_title(
            'Critical $\lambda$ value as a function of the inverse step $\Delta \lambda^{-1}$'
            )

    plt.show()
    try:
        raise exception
    except UnboundLocalError:
        pass

    return inverse_lambda_step_array, lambda_critical_array

    
