def plot_rho_A0(
        self,
        plot_rho,
        ax_rho,
        plot_A0_induced,
        ax_A0_induced,
        alpha,
        ):

    if plot_rho:
        ax_rho.plot(
                self.z,
                self.rho,
                'b',
                alpha=alpha
                )

    if plot_A0_induced:
        ax_A0_induced.plot(
                self.z,
                self.A0_induced,
                'r',
                alpha=alpha
                )
