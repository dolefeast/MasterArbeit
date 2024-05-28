from itertools import count
from numpy import sqrt, mean

def root_mean_square(x):
    return sqrt(mean(x**2))

def update_eigenstates_converge(
        self,
        verbose: int=1, 
        plot_rho: bool=False,
        ax_rho=None,
        plot_A0_induced: bool=False,
        ax_A0_induced=None,
        tol=1e-2,
        max_nodes=5e6,
        ):
    """
    Performs the update_eigenstates until convergence to some tolerance tol
    Parameters:
        verbose: 0,1,2, whether the progress is to be printed.
            if verbose=0, nothing is printed.
            if verbose=1, just the starting message.
            if verbose=2, also the progress at each step.
        plot_rho: bool=False, whether the charge density rho is to be plotted at each iteration.
        ax_rho: None, The maptlotlib.pyplot.ax for the rho to be plotted to.
        plot_A0: bool=False, whether the A0_induced is to be plotted at each iteration.
        ax_A0_induced: None, The maptlotlib.pyplot.ax for A0_induced to be plotted to.
    """
    if verbose:
        print(f"Iterating until convergence")

    # Initializing the calculation
    if verbose == 2:
            print(f'Starting iteration nº 1')

    self.update_eigenstates(
            max_nodes=max_nodes,
            )
    previous_rho = self.rho.copy()
    if self.broken:
        return

    alpha = 0.3
    self.plot_rho_A0(
            plot_rho,
            ax_rho,
            plot_A0_induced,
            ax_A0_induced,
            alpha,
            )

    # To calculate the difference between steps
    # and measuring convergence
    if verbose == 2:
        print(f'Starting iteration nº 2')
    self.update_eigenstates(
            max_nodes=max_nodes,
            )

    alpha = 0.4
    self.plot_rho_A0(
            plot_rho,
            ax_rho,
            plot_A0_induced,
            ax_A0_induced,
            alpha,
            )

    r = root_mean_square(
            (self.rho - previous_rho)
            / ( previous_rho  )
            )

    for iteration in count(3):
        alpha = 0.5
        if r <= tol:
            print(f'Convergence was reached after {iteration-1} iteration{"s"*bool(iteration-1)}')
            break
        previous_rho = self.rho.copy()

        if verbose==2:
            print(f'Starting iteration nº {iteration}')

        self.update_eigenstates(
                max_nodes=max_nodes,
                )
        if self.broken:
            print("The calculation broke")
            break
    
        r = root_mean_square(
                (self.rho - previous_rho)
                / ( previous_rho  )
                )

        self.plot_rho_A0(
            plot_rho,
            ax_rho,
            plot_A0_induced,
            ax_A0_induced,
            alpha,
            )

    # If the iteration stopped, but broken==False,
    # plot the solutoin
    if not self.broken:
        alpha = 1
        self.plot_rho_A0(
            plot_rho,
            ax_rho,
            plot_A0_induced,
            ax_A0_induced,
            alpha,
            )
