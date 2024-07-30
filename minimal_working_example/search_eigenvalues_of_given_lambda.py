from mpmath import pcfd as pbdv
import numpy as np

import scipy.special as special
from scipy.optimize import newton

PI = np.pi

def y1(a,z):
    return np.exp(-0.25*(z**2.0))*special.hyp1f1(0.5*a+0.25,0.5,0.5*(z**2.0))

def y2(a,z):
    return z*np.exp(-0.25*(z**2.0))*special.hyp1f1(0.5*a+0.75,1.5,0.5*(z**2.0))

def U(a,z):
    zeta = 0.5*a+0.25
    return (1/np.sqrt(PI))*(1/(2.0**zeta))*(np.cos(PI*zeta)*special.gamma(0.5-zeta)*y1(a,z) \
    -np.sqrt(2)*np.sin(PI*zeta)*special.gamma(1-zeta)*y2(a,z))

def pbdv(v,z):
    b = -v-0.5
    return U(b,z)

def pbdv_derivative(v, z):
    return z * pbdv(v, z) - pbdv(v+1, z)

def zbar(z, omega, lambda_value):
    return (
            (1 + 1j) / np.sqrt(lambda_value)
            * (omega 
                + lambda_value 
                * (z - 1/2)
                )
            )

def Klein_Gordon_general_solution(
        z,
        omega,
        lambda_value,
        a,
        m=0,
        ):

    if m!=0:
        raise ValueError("m!=0 is not considered")

    nu = -1/2 # Introducing complex values of nu is not easy
              # therefore only study massless case

    first_term = pbdv(
            nu,
            zbar(z, omega, lambda_value),
            )

    second_term = pbdv(
            nu.conjugate(),
            - zbar(z, omega, lambda_value).conjugate(),
            )

    # Watch out! Only applicable for Dirichlet boundary conditions
    b = - a * pbdv(nu, zbar(0, omega, lambda_value)) / pbdv(nu.conjugate(), -zbar(0, omega, lambda_value).conjugate())
    
    return  a * first_term +  b * second_term


def phi_diff_omega(
        omega,
        lambda_value,
        a,
        ):
    # Calculate the derivative from 
    # https://functions.wolfram.com/HypergeometricFunctions/ParabolicCylinderD/20/ShowAll.html
    # with respect to omega, at z = 1
    nu = -1/2

    a_pbdv_at_0 = pbdv(nu, zbar(0, omega, lambda_value))
    b_pbdv_at_0 = pbdv(nu, -zbar(0, omega, lambda_value).conjugate())
    a_pbdv_at_1 = pbdv(nu, zbar(0, omega, lambda_value))
    b_pbdv_at_1 = pbdv(nu, -zbar(0, omega, lambda_value).conjugate())

    b = - a * pbdv(nu, zbar(0, omega, lambda_value)) / pbdv(nu.conjugate(), -zbar(0, omega, lambda_value).conjugate())

    # From the b derivative
    inside_big_parentheses = (
        (1+1j)  # The derivative of z wrt omega
        / np.sqrt(lambda_value)
        *(
            pbdv_derivative(nu, zbar(0, omega, lambda_value))
            ) * pbdv(nu, -zbar(0, omega,lambda_value).conjugate())
        + (1-1j) # The derivative of z wrt omega
        / np.sqrt(lambda_value)
        * (
            pbdv(nu, zbar(0, omega, lambda_value))
            * pbdv_derivative(nu, -zbar(0, omega, lambda_value).conjugate())
        )
        )
    big_derivative_without_a = (
        - a
        / pbdv(nu, -zbar(0, omega, lambda_value)) ** 2
        * inside_big_parentheses
        * pbdv(nu, -zbar(0, omega, lambda_value).conjugate())
        )

    first_term_without_a =  (
        (1+1j)
        / np.sqrt(lambda_value)
        * pbdv_derivative(nu, zbar(1, omega, lambda_value))
        )

    second_term_without_b =  (
        (-1+1j)
        / np.sqrt(lambda_value)
        * pbdv_derivative(nu, -zbar(1, omega, lambda_value).conjugate())
        )
    
    return (
        a * first_term_without_a 
        + b * second_term_without_b
        + a * big_derivative_without_a
        )


def search_eigenvalues(
        omega_guess,
        lambda_value,
        a,
        ):
    #newton(func, x0, fprime=None, args=(), tol=1.48e-08, maxiter=50, fprime2=None, x1=None, rtol=0.0, full_output=False, disp=True)
    # Mission: To find the omega_n for which phi(1, omega)=0

    omega = newton(
            lambda omega, *args: Klein_Gordon_general_solution(1, omega, *args),
            omega_guess,
            args=(lambda_value, a),
            )
    
    return omega

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    n_points = 100

    omega = 1
    lambda_value = 14
    a = 0.01
    b = 1

    z = np.linspace(0, 1, n_points)

    phi = lambda z: Klein_Gordon_general_solution(
            z,
        omega,
        lambda_value,
        a,
        b,
        m=0,
        )

    omega = search_eigenvalues(np.pi, lambda_value=lambda_value, a=a)
    print(omega)


    # plt.plot(z, phi(z))
    plt.show()


