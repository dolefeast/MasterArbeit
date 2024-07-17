from cycler import cycler

def plot_rho_A0(
        self,
        plot,
        ax_rho,
        color_rho,
        ax_A0_induced,
        color_A0_induced,
        alpha,
        ):

    if plot:
        ax_rho.plot(
                self.z,
                self.rho,
                alpha=alpha,
                color=color_rho,
                )
        ax_A0_induced.plot(
                self.z,
                self.A0_induced,
                alpha=alpha,
                color=color_A0_induced,
                )
