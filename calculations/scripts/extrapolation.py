def linear_extrapolation(x:[float], y:[float], x_step=None):
    """Returns the next step of an array y by linear extrapolation (Euler method)"""

    assert sorted(x) == x
    dy = (y[-1] - y[-2])
    dx = (x[-1] - x[-2])

    if x_step is None:
        x_step = dx
    return x[-1] + x_step, y[-1] + dy/dx * x_step
