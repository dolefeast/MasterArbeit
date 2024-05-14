def update_eigenstates_script(
        self,
        n_iterations:int=None,
        save_solutions: bool=False,
        verbose: int=1,
        tol:float=1e-2,
        plot_rho: bool=False,
        ax_rho=None,
        plot_A0_induced: bool=False,
        ax_A0_induced=None,
        max_nodes=5e6,
        ):
    """
    Parameters:
        n_iterations:int=None, if None it iterates until convergence is achieved. If an int,
            it iterates a fixed amount of times.
        save_solutions: bool=False, whether if the solutions are to be saved
        verbose: int=1, can be 0, 1, 2. Controls the verbosity of the scripts
        tol=float:1e-2, the tolerance to which convergence is considered
    """

    if n_iterations is None:
        self.update_eigenstates_converge(
                verbose=verbose,
                plot_rho=plot_rho,
                ax_rho=ax_rho,
                plot_A0_induced=plot_A0_induced,
                ax_A0_induced=ax_A0_induced,
                tol=tol,
                max_nodes=max_nodes,
                )

    elif isinstance(n_iterations, int):
        self.update_eigenstates_iterate(
                n_iterations=n_iterations,
                verbose=verbose,
                plot_rho=plot_rho,
                ax_rho=ax_rho,
                plot_A0_induced=plot_A0_induced,
                ax_A0_induced=ax_A0_induced,
                max_nodes=max_nodes,
                )

    if save_solutions and not self.broken:
        self.save_solutions()
