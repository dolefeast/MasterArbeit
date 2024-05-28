def main(
        m=0,
        directory="",
        a_values=[],
        ):
    from plot.format_omega_evolution import format_omega_evolution
    from plot.plot_eigenvalue_evolution import plot_eigenvalue_evolution

    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    if isinstance(a_values, (int, float)):
        a_values = [a_values]
    elif len(a_values)==0:
        a_values = range(10)

    format_omega_evolution(ax)

    for a in a_values:
        plot_eigenvalue_evolution(
                m,
                a,
                ax=ax,
                directory=directory,
                )

    plt.show()

