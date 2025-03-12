import numpy as np
import os
from itertools import cycle
from scipy.interpolate import CubicSpline
from pathlib import Path
from mpmath import quad
import shutil

class Vacuum_Polarization:
    def __init__(self,
            maxN=1,
            m=0,
            a=1,
            e=1,
            bcs="dirichlet",
            lambdaMin=10,
            lambdaMax=25,
            lambdaStep=2,
            nIterations=None,
            pointsPerMode = 8
            ):

        self.e = e
        self.a = a
        self.maxN = maxN
        self.m = m
        self.lambdaValue = lambdaMin
        self.bcs = bcs


        if maxN<5:
            self.nPoints = 200
            self.ambjorn = True # Controls whether the self.e**2/np.pi * A0(z) term is added at the end
            self.smoothing = False
        else:
            print("Watch out, the number of points is modified!")
            if pointsPerMode != 8:
                print("Warning: To get a proper filtering, pointsPerMode should be 8")

            self.nPoints = pointsPerMode * ( maxN + 1 )
            self.ambjorn = ambjorn      # This can be changed depending on the results we are looking for
            self.smoothing = smoothing  # This can be changed depending on the results we are looking for
        
        # The mesh points
        self.z = np.linspace(0, 1, self.nPoints)
        
        self.lambdaMin = lambdaMin
        self.lambdaMax = lambdaMax
        self.lambdaStep = lambdaStep
        self.lambdaArray = []

        self.rho = [0] * self.nPoints
        self.A0Induced = lambda z: 0
        self.A0 = lambda z: - self.lambdaValue * (z - 1/2) 

    from Vacuum_Polarization.physics.KleinGordonEquation import KleinGordonEquation
    from Vacuum_Polarization.physics.perturbativeApproximations import (
            perturbativePhi,
            perturbativeModeRho,
            perturbativeTotalVacuumPolarization
            )

    from Vacuum_Polarization.numerics.bisectionMethod import (
            bisectionMethod,
            eigenvaluesUpperBound,
            eigenvaluesLowerBound
            )
