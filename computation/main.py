from snippets.eigenvalue_evolution_wrt_a import main

main(
    m=0,
    E_count=30,
    a_min=0.001,
    a_max=0,5,
    a_count=5,
    n_iterations=None,
    verbose=0,
    tol=5e-6,
    max_nodes=5e5,
    plot=False,
    save_solutions=True,
    read_solutions=False,
        )

