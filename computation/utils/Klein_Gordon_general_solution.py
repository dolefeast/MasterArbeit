from mpmath import pcfd as pbdv
import numpy as np

import numpy
import scipy.special as special

PI = numpy.pi

def y1(a,z):
    return numpy.exp(-0.25*(z**2.0))*special.hyp1f1(0.5*a+0.25,0.5,0.5*(z**2.0))

def y2(a,z):
    return z*numpy.exp(-0.25*(z**2.0))*special.hyp1f1(0.5*a+0.75,1.5,0.5*(z**2.0))

def U(a,z):
    zeta = 0.5*a+0.25
    return (1/numpy.sqrt(PI))*(1/(2.0**zeta))*(numpy.cos(PI*zeta)*special.gamma(0.5-zeta)*y1(a,z) \
    -numpy.sqrt(2)*numpy.sin(PI*zeta)*special.gamma(1-zeta)*y2(a,z))

def pbdv(v,z):
    b = -v-0.5
    return U(b,z)

def Klein_Gordon_general_solution(
        z,
        lambda_value,
        omega,
        a,
        b,
        m=0,
        ):

    nu = -1/2 # Only worrying right now of m = 0
    zbar = (
            (1 + 1j) / np.sqrt(lambda_value)
            * (omega 
                + lambda_value 
                * (z - 1/2)
                )
            )

    first_term = pbdv(
            nu,
            zbar,
            )

    second_term = pbdv(
            nu.conjugate(),
            - zbar.conjugate()
            )

    b = - a * first_term[0] / second_term[0]
    
    return  a * first_term +  b * second_term

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    z  = np.linspace(0, 1, 400)
    lambda_value = 16.5
    omega = 0.
    a = 2
    b = 1

    plt.plot(
            z,
            Klein_Gordon_general_solution(
        z=z,
        lambda_value=lambda_value,
        omega=omega,
        a=a,
        b=b,
        )
        )

    plt.show()
