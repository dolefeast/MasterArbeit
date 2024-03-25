from math import pi, sin, cos, sqrt
from scripts.iterable_output import iterable_output

def sign(x):
    if x>0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0

@iterable_output
def dirichlet_eigenstate(z, omega_n, m, lambda_value):
    n = sign(omega_n)*sqrt(abs(omega_n**2 - m**2))/pi

    if n == 0:
        return 0

    eigenstate= (
            abs(omega_n) ** -(1/2)
            *(
                sin(pi * n * z) 
                + lambda_value  *
                omega_n 
                / (2 * pi * abs(n))
                * ((1/2-z)/pi/n * sin(pi * n * z)
                    - z * (1-z) * cos(pi * n * z)))
                )
        
    return eigenstate

@iterable_output
def dirichlet_eigenstate_gradient(z, omega_n, m, lambda_value):
    n = sign(omega_n)*sqrt(abs(omega_n**2 - m**2))/pi

    if n == 0:
        return 0

    eigenstate= (
            abs(omega_n) ** -(1/2)
            *(
                cos(pi * n * z) * pi*n
                + lambda_value  *
                omega_n *(
                 (-1/pi/n * sin(pi * n * z)
                +(1/2-z)*cos(pi * n * z) # Deriv of first term
                    - (1-z) * cos(pi * n * z)
                    + z * cos(pi * n * z))
                    + pi*n*z * (1-z) * sin(pi * n * z))
                / (2 * pi * abs(n))
                )
            )
        
    return eigenstate


if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    m = 5
    lambda_value = 1
    z = np.linspace(0, 1, 100)

    n = np.arange(-10, 10, 1)
    omega_array = np.sign(n) * np.sqrt(n**2 * np.pi ** 2 + m**2)

    fig, (ax1, ax2) = plt.subplots(1, 2)
    
    for omega_n in omega_array:
        ax1.plot(z, dirichlet_eigenstate(z, omega_n, m, lambda_value), 'b-')
        ax2.plot(z, dirichlet_eigenstate_gradient(z, omega_n, m, lambda_value), 'g-')
        ax2.plot(z, dirichlet_eigenstate_gradient(z, omega_n, m, lambda_value), 'g-')

    plt.show()

        
