def magnitude_evolution(
        m,
        a,
        ax,
        magnitude='rho',
        bcs='dirichlet',
        lambda_value=None,
        sig_digs=3,
        color='b',
        directory="",
        max_alpha=1,
        ):
    """
    Parameters,
        m: float, the mass corresponding to the magnitude we want to plot
        a: float, the interval size corresponding to the magnitude we want to plot
        magnitude: float, the magnitude we want to plot
        bcs: str='dirichlet', the boundary conditions used to calculate said magnitude
        lambda_value, may be scalar or list of scalars. The field strength of the magnitude we want
        sig_digs, the significative digits to which to read the lambda_value

    Adds plots to an existing axis.
    """

    from plot.plot_single_magnitude import plot_single_magnitude
    from utils.read_files_fixed_m_a import read_files_fixed_m_a
    from utils.get_lambda_index import get_lambda_index
    from utils.plot_from_0_to_1 import plot_from_0_to_1

    data = read_files_fixed_m_a(
            m=m,
            a=a,
            read_things=magnitude,
            bcs='dirichlet',
            sig_digs=sig_digs,
            directory=directory,
            )

    # If the given lambda_value is not in the calculated lambda values array
    # then search for the closest one available that is smaller than the given one
    exact_lambda_value = False
    factor = 1
    if magnitude ==  "A0_induced" and m == 0:
        factor = 1/a
    try:
        # Check if lambda_value is iterable
        lambda_value_array = list(lambda_value)
    except TypeError:
        if isinstance(lambda_value, (int, float)):
            lambda_value_array = [lambda_value]
        elif lambda_value is None:
            lambda_value_array = data["lambda_value"]
            # These are calculated
            exact_lambda_value = True

    max_lambda = max(lambda_value_array)

    # I need the index to be able to plot the desired magnitude
    for i, lambda_value in enumerate(lambda_value_array):
        if not exact_lambda_value:
            i = get_lambda_index(data["lambda_value"], lambda_value)

        alpha = (
                0.2 + 0.8 * (lambda_value / max_lambda) ** 4
                ) * max_alpha
        ax.plot(*plot_from_0_to_1(data[magnitude][i] * factor),
                color,
                alpha=alpha,
                )
