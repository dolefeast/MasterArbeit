def main(
        m,
        lambda_min,
        n_iterations,
        verbose,
        save_solutions,
        read_solutions,
        max_nodes=5e8,
        inverse_lambda_step_min = 5,
        inverse_lambda_step_max = 25,
        inverse_lambda_step_count = 10,
        tol=1e-4,
        ):

    import matplotlib.pyplot as plt
    from numpy import linspace

    from physics import Vacuum_Polarization
    from tests.create_new_class import main as create_new_class_main

    fig, ax = plt.subplots()

    inverse_lambda_step_array = linspace(
        inverse_lambda_step_min,
        inverse_lambda_step_max,
        inverse_lambda_step_count,
            )

    lambda_critical_array = []

    try:
        for inverse_lambda_step in inverse_lambda_step_array:
            lambda_critical = create_new_class_main(
                m,
                lambda_min,
                inverse_lambda_step**-1,
                n_iterations,
                verbose,
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

    
