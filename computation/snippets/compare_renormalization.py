import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

from plot.compare_renormalization import compare_renormalization
from utils.read_files_fixed_m_a import read_files_fixed_m_a

def main(
        lambda_value: float,
        m=0,
        a=1,
        bcs='dirichlet',
        sig_digs=3,
        ax_rho=None,
        ax_A0_induced=None,
        directory="",
        max_alpha=1,
        ):

    # If nothing is said about both axes, generate them
    show = False
    if ax_rho is None:
        import matplotlib.pyplot as plt
        fig, ax_rho = plt.subplots(1)

        ax_rho.set_ylabel(r"$\rho(z)$")
        ax_rho.set_xlabel(r"$z$")
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


    rho, mode_sum_rho = compare_renormalization(
            a,
            lambda_value,
            data,
            ax_rho,
            )

    if show:
        plt.show()
