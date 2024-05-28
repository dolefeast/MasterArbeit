def main(
        m=0,
        a=1,
        lambda_value=None,
        bcs='dirichlet',
        sig_digs=3,
        ax_rho=None,
        ax_A0_induced=None,
        directory="",
        max_alpha=1,
        ):

    from plot.magnitude_evolution import magnitude_evolution

    # If nothing is said about both axes, generate them
    show = False
    if ax_rho is None and ax_A0_induced is None:
        import matplotlib.pyplot as plt
        fig, (ax_rho, ax_A0_induced) = plt.subplots(2)

        ax_rho.set_ylabel(r"$\rho(z)$")
        ax_rho.set_xlabel(r"$z$")
        ax_A0_induced.set_ylabel(r"$\tilde{A}_0(z)$")
        ax_A0_induced.set_xlabel(r"$z$")
        show = True

    if isinstance(a, (int, float)):
        magnitude_evolution(
                m, 
                a, 
                ax_rho,
                magnitude='rho',
                lambda_value=lambda_value,
                color='b',
                directory=directory,
                max_alpha=max_alpha,
                )

        magnitude_evolution(
                m, 
                a, 
                ax_A0_induced,
                magnitude='A0_induced',
                lambda_value=lambda_value,
                color='r',
                directory=directory,
                max_alpha=max_alpha,
                )
    else:
        for a_value in a:
            try:
                main(
                    m=m,
                    a=a_value,
                    lambda_value=lambda_value,
                    bcs=bcs,
                    sig_digs=sig_digs,
                    ax_rho=ax_rho,
                    ax_A0_induced=ax_A0_induced,
                    directory=directory,
                    max_alpha=max_alpha,
                )
            except ValueError:
                print(f"a={a} was empty")
                continue


    if show:
        plt.show()
