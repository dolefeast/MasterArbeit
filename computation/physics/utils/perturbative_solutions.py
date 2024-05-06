from math import pi, sin, cos, sqrt
from utils.iterable_output import iterable_output
from numpy import sign

@iterable_output
def dirichlet_eigenstate(z, omega_n, m, lambda_value):
    n = sign(omega_n) * sqrt(abs(omega_n**2 - m**2)) / pi

    if n == 0:
        return 0

    eigenstate = abs(omega_n) ** -(1 / 2) * (
        sin(pi * n * z)
        + lambda_value
        * omega_n
        / (2 * pi * abs(n))
        * ((1 / 2 - z) / pi / n * sin(pi * n * z) - z * (1 - z) * cos(pi * n * z))
    )

    return eigenstate

@iterable_output
def dirichlet_eigenstate_gradient(z, omega_n, m, lambda_value):
    n = sign(omega_n) * sqrt(abs(omega_n**2 - m**2)) / pi

    if n == 0:
        return 0

    eigenstate = abs(omega_n) ** -(1 / 2) * (
        cos(pi * n * z) * pi * n
        + lambda_value
        * omega_n
        * (
            (
                -1 / pi / n * sin(pi * n * z)
                + (1 / 2 - z) * cos(pi * n * z)  # Deriv of first term
                - (1 - z) * cos(pi * n * z)
                + z * cos(pi * n * z)
            )
            + pi * n * z * (1 - z) * sin(pi * n * z)
        )
        / (2 * pi * abs(n))
    )

    return eigenstate


@iterable_output
def neumann_eigenstate(z, omega_n, m, lambda_value):
    n = sign(omega_n) * sqrt(abs(omega_n**2 - m**2)) / pi

    if n == 0:
        return 0

    eigenstate = abs(omega_n)**(-1/2) * (
            cos(pi * n * z)
            + lambda_value 
                * omega_n
                / (2*pi * abs(n))
                * (
                    ( 1/2 - z )
                    * cos(pi * n * z)
                    / pi / n
                    + (
                        z * (1-z) 
                        + ( pi * n ) ** -2
                    * sin( pi * n * z )
                        )
                    )
                )

    return eigenstate

@iterable_output
def neumann_eigenstate_gradient(z, omega_n, m, lambda_value):
    n = sign(omega_n) * sqrt(abs(omega_n**2 - m**2)) / pi

    if n == 0:
        return 0

    eigenstate = abs(omega_n)**(-1/2) * (
            - pi * n * sin(pi * n * z)
            + lambda_value 
                * omega_n
                / (2*pi * abs(n))
                * (
                    n * pi 
                    * z * (1-z) 
                    * cos( n * pi * z)
                    + (1/2 - z) * sin(n * pi * z)
                    )
                )

    return eigenstate
