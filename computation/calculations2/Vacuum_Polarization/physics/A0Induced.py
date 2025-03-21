from scipy.interpolate import CubicSpline

def calculateA0Induced(self):
    rhoInterpolated = CubicSpline(self.z, self.rho)

    A0InducedShifted = rhoInterpolated.antiderivative(2)
    offset = A0InducedShifted(1/2)
    self.A0Induced = lambda z: -(A0InducedShifted(z) - offset)
