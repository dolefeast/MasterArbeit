from tests.create_new_class import main

main(
    m = 0,
    lambda_min = 0.5,
    lambda_step = 1,
    n_iterations = 3,
    verbose = 1,
    tol = 5e-8,
    max_nodes=5e5,
    plot=True,
    # plot_rho=True,
    # plot_A0_induced=True,
    save_solutions=False,
    read_solutions=True,
    coupling = 5,
        )

