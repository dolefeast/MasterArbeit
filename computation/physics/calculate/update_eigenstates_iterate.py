def update_eigenstates_iterate(
        self,
        n_iterations: int,
        verbose: int=1, 
        plot_rho: bool=False,
        ax_rho=None,
        plot_A0_induced: bool=False,
        ax_A0_induced=None,
        max_nodes=5e6,
        ):
    """
    Performs the update_eigenstates method n_iterations times
    Parameters:
        n_iterations: int, the number of iterations to be performed
        verbose: 0,1,2, whether the progress is to be printed
            if verbose=0, nothing is printed
            if verbose=1, just the starting message
            if verbose=2, also the progress at each step
        plot_rho: bool=False, whether the charge density rho is to be plotted at each iteration
        ax_rho: None, The maptlotlib.pyplot.ax for the rho to be plotted to
        plot_A0: bool=False, whether the A0_induced is to be plotted at each iteration
        ax_A0_induced: None, The maptlotlib.pyplot.ax for A0_induced to be plotted to
    """
    if verbose:
        print(f"Calculating {n_iterations} iteration{'s'*bool(n_iterations-1)}")

    for iteration in range(n_iterations):
        if verbose==2:
            print(f'Starting iteration nº {iteration+1}')
        self.update_eigenstates(
                max_nodes=max_nodes,
                )

        if self.broken:
            print("The calculation broke")
            break

        # To distinguish between different solutions
        alpha = (
                0.3 
                + 0.7*(
                    (iteration+1) / n_iterations 
                    ) **3
                )

        self.plot_rho_A0(
                plot_rho,
                ax_rho,
                plot_A0_induced,
                ax_A0_induced,
                alpha,
                )
