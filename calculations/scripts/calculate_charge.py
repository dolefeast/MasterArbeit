import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad

def calculate_charge(rho, kind='linear'):
    """Calculates the induced charge from 0 to 1/2"""
    z = np.linspace(0, 1, len(rho))
    rho = interp1d(z, rho, kind=kind)

    rho_callable = lambda z: rho(z)

#    q = mpmath.quad(rho_callable, [0,1/2])
    q = quad(rho_callable, 0,1/2)[0]

    return q

if __name__ == "__main__":
    x = np.linspace(0, 1, 100)
    rho = -0.05 * np.sin(2*np.pi * x)


    print(calculate_charge(rho)[0])

