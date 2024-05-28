import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

from plot.compare_renormalization import compare_renormalization
from utils.read_files_fixed_m_a import read_files_fixed_m_a
from utils.plot_from_0_to_1 import plot_from_0_to_1

def main(
        lambda_value: float,
        m=0,
        a=1,
        bcs='dirichlet',
        sig_digs=3,
        ax_rho=None,
        ax_E_induced=None,
        ax_A0_induced=None,
        directory="",
        max_alpha=1,
        ):

    # If nothing is said about both axes, generate them
    show = False
    if ax_rho is None and ax_E_induced is None and ax_A0_induced is None:
        import matplotlib.pyplot as plt
        fig1, ax_rho = plt.subplots(1)
        fig2, ax_E_induced = plt.subplots(1)
        fig3, ax_A0_induced = plt.subplots(1)

        ax_rho.set_ylabel(r"$\rho(z)$")
        ax_E_induced.set_ylabel(r"$\tilde{E}$")
        ax_A0_induced.set_ylabel(r"$\tilde{A}_0$")
        ax_rho.set_xlabel(r"$z$")
        ax_E_induced.set_xlabel(r"$z$")
        ax_A0_induced.set_xlabel(r"$z$")
        show = True

    data = read_files_fixed_m_a(
            m,
            a,
            read_things=[
                'rho',
                'A0_induced',
                ],
            bcs=bcs,
            sig_digs=sig_digs,
            directory=directory,
            )


    z, rho, rho_mode_sum = compare_renormalization(
            a,
            lambda_value,
            data,
            ax_rho,
            )

    rho = CubicSpline(
            z,
            rho,
            )

    rho_mode_sum = CubicSpline(
            z,
            rho_mode_sum,
            )

    E_induced = rho.antiderivative(1)
    E_induced_mode_sum = rho_mode_sum.antiderivative(1)

    ax_E_induced.plot(z, E_induced(z))
    ax_E_induced.plot(z, E_induced_mode_sum(z))

    if show:
        fig1.tight_layout()
        fig2.tight_layout()
        fig3.tight_layout()
        plt.show()
