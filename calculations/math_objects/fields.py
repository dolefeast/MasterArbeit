"""
The ambient will be governed by two fields, the classical
electric field and the scalar field. 
"""
from __init__ import *
from dataclasses import dataclass
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np
import sys


@dataclass
class Field(object):
    """
    The properties of the scalar field in question:
    Parameters:
        mass: The mass in natural units
        charge: Charge in mass units (since 1+1 dimensions)

    Properties:
        z: the space [0, 1] considered
        field: the field at each point z
        deriv: the derivative wrt z of the field
    """

    # Implement here Robin boundary conditions?
    value: None = None
    mass: float = 0
    charge: float = 0
    n_points: int = 50
    _gradient: np.ndarray = None
    _conjugate: np.ndarray = None
    name: str = "Field"

    def __post_init__(self):
        if self.value is None:
            self.value = np.zeros(self.n_points)  # If value is not given, create
            # empty sequence.
        else:
            try:
                iter(self.value)
                self.value = np.array(self.value)  # If value is an array, leave as is.
                # Change into numpy array just in
                # case.
                self.n_points = len(self.value)
            except TypeError:
                self.value = self.value * np.ones(
                    self.n_points
                )  # If value is scalar, it
                # was expected to be constant
                # value, thus create new array
        self.z = np.linspace(0, 1, self.n_points)
        self.callable = None

    def __repr__(self):
        return f"{self.name}(mass={self.mass}, charge={self.charge})"

    def assertion(self, other_field):
        assert self.mass == other_field.mass, "The masses of the fields do not match"
        assert (
            self.charge == other_field.charge
        ), "The masses of the fields do not match"
        assert (
            self.n_points == other_field.n_points
        ), "n_points of the fields do not match"

    def __add__(self, other_field):
        self.assertion(other_field)

        return Field(
            mass=self.mass,
            charge=self.charge,
            n_points=self.n_points,
            value=self.value + other_field.value,
        )

    def __mul__(self, k):
        try:
            iter(k)
            is_iterable = True
        except TypeError:
            is_iterable = False
        except AssertionError:
            raise AssertionError("The shapes of the scalar and the field do not match")
        if isinstance(k, (float, int)) or isinstance(k, complex) or is_iterable:
            return Field(
                mass=self.mass,
                charge=self.charge,
                n_points=self.n_points,
                value=k * self.value,
            )
        else:
            raise TypeError("k must be a real or complex scalar (function)!")
        pass

    def __rmul__(self, k):
        return self * k  # Multiplication is in this case commutative

    def __call__(self, z):
        """returns field(z) with z between 0 and 1 by interpolation."""
        try:
            z_min, z_max = min(z), max(z)  # means z was an iterable
            assert z_min >= 0, f"z_min={z_min}. Every item in z must be greater than 0."
            assert z_max <= 1, f"z_max={z_max}. Every item in z must be smaller than 1."
        except TypeError:
            assert (
                z >= 0 and z <= 1
            ), f"z={z}. z must be greater than 0 and smaller than 1."

        if self.callable is None:
            self.callable = self.calculate_interpolation()
        return self.callable(
            z
        )  # callable accepting arrays is already built in the PPoly object

    def calculate_interpolation(self):
        """Calculates interpoltion of the field value property"""
        field_interpolation = sp.interpolate.CubicSpline(self.z, self.value)
        return field_interpolation

    def calculate_field_integral(self, method="quad"):
        """Calculates the integral of the field between 0 and 1. Admitted scipy integration methods, 'quad', 'fixed_quad', 'guassian', 'quadrature'."""
        methods = {
            "quad": sp.integrate.quad,
            "fixed_quad": sp.integrate.fixed_quad,
            "gaussian": sp.integrate.fixed_quad,  # As from integrate docs, it is repeated
            "quadrature": sp.integrate.quadrature,
        }

        z0 = min(self.z)  # Should be 0
        z1 = max(self.z)  # Should be 1
        if method not in methods:
            raise TypeError(
                f"The supplied method of integration {method} is not a valid one!\n\t Please input one of the following: {self.methods.keys()}"
            )

        return self.methods[method](self, z0, z1)[
            0
        ]  # Watch out! This is total integral of the field, not F(x) s.t. F'(x) = phi(x)

    def norm_squared(self, method="quad"):
        """Calculates the square of the L2 norm of the field
        Parameters:
            method: str. The method used for the integration in scipy.integrate
        Returns:
            norm_squared: float. The norm squared of the field. The integral from 0 to 1 of the wave function times itself, complex conjugated
        """
        return sp.integrate.quad(self.value * np.conjugate(self.value))

    def normalize(self, method="quad"):
        norm = np.sqrt(self.norm_squared)
        self.value /= norm

    def conjugate(self):
        return np.conjugate(self.value)

    @property
    def gradient(self):  # To calculate the derivative each time
        # the object is generated is costly
        # and absurd. This way it will only
        # generate it if called.
        if self._gradient is None:
            self._gradient = self.calculate_gradient()

        return self._gradient

    def calculate_gradient(self, order=1):
        diff_phi = np.diff(self.value, order)
        diff_phi = np.append(diff_phi, diff_phi[-1])
        # Due to chain rule we are missing a factor n_point
        # diff_phi = diff_phi * self.n_points
        Diff_phi = Field(
            value=diff_phi,
            mass=self.mass,
            charge=self.charge,
            n_points=self.n_points,
            name=self.name + " gradient",
        )
        return Diff_phi


class Vector_Potential(Field):
    """The electric field assuming static solutions comes from
    A_mu = (A_0, 0), with -nablaA_0 = E. Using Gauß gauge,
    for an initial constant electric field we get
    A0 = -lambda (z - 1/2)"""

    # As of meeting 15.11.23, note that the electric field
    # will be generated by some charges outside of the
    # space we are considering. This means that independently
    # of whatever happens inside the region z \in [0, 1], the
    # electric field MUST be continuous and thus equal to
    # some fixed value on the boundaries.

    _electric_field: np.ndarray = None

    @property
    def electric_field(self):  # To calculate the derivative each time
        # the object is generated is costly
        # and absurd. This way it will only
        # generate it if called.
        if self._gradient is None:
            self._gradient = self.calculate_gradient()

        return self._gradient

    def calculate_gradient(self, order=1):
        diff_phi = np.diff(self.value, order)
        diff_phi = np.append(diff_phi, diff_phi[-1])
        return diff_phi * self.n_points


if __name__ == "__main__":
    n_points = 500

    # f = lambda z, sigma, z0: 1/np.sqrt(2*np.pi*sigma**2) * np.exp(-(z-z0)**2/2/sigma**2)
    f = lambda z, sigma, z0: z**2
    g = lambda z, sigma, z0: np.cos(np.pi * z)
    parabola1 = f(np.linspace(0, 1, n_points), 0.025, 0.5)
    parabola2 = g(np.linspace(0, 1, n_points), 0.025, 0.5)
    phi1 = Field(value=parabola1, name="phi1")
    phi2 = Field(value=parabola2, name="phi2")

    plt.plot(phi2.z, phi2.value, label=phi2)
    plt.plot(phi2.z, phi2.gradient(phi2.z), label=phi2)

    plt.legend(loc="best")

    plt.show()
