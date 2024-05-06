import numpy as np
from scipy.interpolate import CubicSpline

class Field:
    def __init__(
            self,
            m: float=0,
            e: float=0,
            value: [float]=0,
            n_points: int=None,
            ):
        self.e = e
        self.m = m

        self.callable = None



        if not isinstance(n_points, (int)):
            raise ValueError("n_points was not int")
        else:
            self.n_points = n_points
            self.z = np.linspace(
                    0,
                    1, 
                    n_points,
                    )

        try:
            # Check if the desired value of the Field is iterable
            float(value)
            self.value = np.array([value] * self.n_points)
        except TypeError:
            self.value = np.array(value)
            assert len(self.value) == len(self.z)

    def __repr__(self):
        return f"{self.value}"

    def __add__(self, other_field):
        return Field(
                m=self.m,
                e=self.e,
                n_points=self.n_points,
                value=self.value + other_field.value
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
            m=self.m,
            e=self.e,
            n_points=self.n_points,
            value=k * self.value,
            )
        else:
            raise TypeError("k must be a real or complex scalar (function)!")
        pass

    def __rmul__(self, k):
        return self * k  # I want multiplication to be commutative

    # To make the object callable
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

        # callable was never calculated
        if self.callable is None:
            self.callable = self.calculate_interpolation()
        # it was, return the desired value
        return self.callable(
            z
            )  # callable accepting arrays is already built in the PPoly object

    def calculate_interpolation(self):
        """Calculates interpoltion of the field value property"""
        field_interpolation = CubicSpline(self.z, self.value)
        return field_interpolation
