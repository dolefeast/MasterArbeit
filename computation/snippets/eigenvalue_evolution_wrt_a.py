def main(
    m=0,
    E_min=0.5,
    E_count=40,
    a_min=0.01,
    a_max=1,
    a_count=10,
    n_iterations=3,
    verbose=1,
    tol=5e-8,
    max_nodes=5e5,
    plot=True,
    # plot_rho=True,
    # plot_A0_induced=True,
    save_solutions=False,
    read_solutions=True,
        ):
    import numpy as np

    from physics import Vacuum_Polarization

    system=Vacuum_Polarization(
            m=0,
            E=E_min,
            read_solutions=False,
            )

    system.save_solutions_dir = "saved_solutions/a_evolution"

    a_array = np.linspace(
            a_min,
            a_max,
            a_count
            )

    # The limit without backreaction used to be (for m=0)
    # lambda = a²E ~ 16
    E_min_array = E_min/a_array**2
    E_max_array = 25/a_array**2

    E_array_array = [
            np.linspace(
                E_min,
                E_max,
                )
            for E_min, E_max in zip(E_min_array, E_max_array)
            ]

    for a, E_array in zip(a_array, E_array_array):
        system.broken = False
        system.a = a
        for E in E_array:
            system.E = E
            system.update_eigenstates_script(
                    n_iterations=n_iterations,
                    verbose=verbose,
                    save_solutions=save_solutions,
                    tol=tol,
                    max_nodes=max_nodes,
                )
