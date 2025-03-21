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
            lambdaStep=1,
            lambdaStepMin=1e-8,
            relaxParameter=1,
            walkbackFactor=0.4,
            directory="",
            saveData=False,
            showPlots=False,
            savePlots=False,
            intermediatePlots=False,
            hadamard=True,
            smoothing=False,
            nIterations=None,
            maxNIterations=100,
            pointsPerMode=8,
            subtractPertModes=False,
            bisectionMaxiter=1500,
            bisectionTol=1e-15,
            convergenceTol=1e-3,
            sigDigs=9,
            parallelization=True,
            extrapolate=True,
            ):

        self.e = e
        self.a = a
        self.maxN = maxN
        self.m = m
        self.lambdaValue = lambdaMin
        self.bcs = bcs

        start = (self.bcs == "dirichlet") # because the neumann boundary conditions admit a n=0 case.
        self.eigenvalues = [ np.sqrt( self.m**2 + (np.pi * n) ** 2 ) for n in range(start, maxN+1) ] 

        if maxN<5:
            self.nPoints = 200
            self.hadamard = False # Controls whether the self.e**2/np.pi * A0(z) term is added at the end
            self.smoothing = False
        else:
            self.hadamard = hadamard      # This can be changed depending on the results we are looking for
            self.smoothing = smoothing  # This can be changed depending on the results we are looking for
            if self.smoothing == True and pointsPerMode != 8:
                raise ValueError("pointsPerMode != 8 and smoothing = True are incompatible settings")
                pointsPerMode = 8


            self.nPoints = pointsPerMode * ( maxN + 1 )
        
        # The mesh points
        self.z = np.linspace(0, 1, self.nPoints)
        
        self.lambdaMin = lambdaMin
        self.lambdaMax = lambdaMax
        self.lambdaStep = lambdaStep
        self.walkbackFactor = walkbackFactor
        self.relaxParameter = relaxParameter
        self.lambdaArray = []

        self.nIterations = nIterations  # = 1 is the external field approximation. = None is the self-consistent solution
        self.maxNIterations = maxNIterations # For safety
        self.lambdaStepMin = lambdaStepMin # For safety

        self.rho = [0] * self.nPoints
        self.A0Induced = lambda z: 0
        self.A0 = lambda z: - self.lambdaValue * (z - 1/2) 

        # Avoids recursion errors
        self.A0InducedHistory = [[self.A0Induced]]
        self.A0History = [[self.A0]]

        self.bisectionMaxiter = bisectionMaxiter
        self.bisectionTol = bisectionTol
        self.convergenceTol = convergenceTol

        self.subtractPertModes = subtractPertModes

        self.saveData = saveData # Stores results as txt
        if saveData and directory == "":
            raise ValueError("If saveData==True, directory cannot be empty")
        self.directory = directory

        if saveData:
            directoryExists = Path("data/"+directory+"/"+bcs).is_dir()
            if directoryExists:
                print(f"Warning: data/{directory} already existed.")
            else:
                Path("data/"+directory+"/"+bcs).mkdir(parents=True, exist_ok=True)
        self.sigDigs = sigDigs

        self.parallelization = parallelization # For the mode calculations
        self.extrapolate = extrapolate  # To predict the folowing A0 after a converged one

        self.colorCycle = cycle(
    "#C52E19FF, #AC9765FF, #54D8B1FF, #B67C3BFF, #175149FF, #AF4E24FF".split(", ")
        )
        self.color = next(self.colorCycle)

        self.showPlots = showPlots
        self.savePlots = savePlots

        self.plotFinal = showPlots or savePlots # Plots converged quantities
        self.intermediatePlots = intermediatePlots and savePlots # These are only saved. Otherwise there'd be too many

        if self.plotFinal:
            self.plotSetup()

        if self.intermediatePlots:
            self.plotIntermediateStepsSetup()

    from Vacuum_Polarization.physics.KleinGordonEquation import KleinGordonEquation
    from Vacuum_Polarization.physics.perturbativeApproximations import (
            perturbativePhi,
            perturbativeModeRho,
            perturbativeTotalVacuumPolarization
            )

    from Vacuum_Polarization.numerics.bisectionMethod import (
            bisectionMethod,
            eigenvaluesUpperBound,
            eigenvaluesLowerBound,
            NoRootFoundError,
            MaxIterError,
            )

    from Vacuum_Polarization.numerics.convergence import (
            initializeA0Lists,
            convergence,
            )

    from Vacuum_Polarization.numerics.singleIteration import (
            singleIteration,
            )

    from Vacuum_Polarization.numerics.fullScript import (
            fullScript,
            walkback
            )

    from Vacuum_Polarization.numerics.calculateEigenstates import (
            calculateNormSquared,
            calculateSingleEigenstate,
            calculateEigenstatesParallel
            )

    from Vacuum_Polarization.physics.A0Induced import (
            calculateA0Induced
            )

    from Vacuum_Polarization.physics.rho import (
            calculateRho,
            convolveRho,
            )

    from Vacuum_Polarization.scripts.saveData import saveDataScript

    from Vacuum_Polarization.scripts.utils import (
            floatToStr,
            strToFloat
            )

    from Vacuum_Polarization.scripts.readData import (
            openSolutionFamilyArray,
            setConfigFromDict,
            readSolutionsFromFile,
            getPosixForQuantities,
            getLambdaValue
            )

    from Vacuum_Polarization.plotting.plotSetup import (
            plotSetup,
            plotIntermediateStepsSetup,
            plotIntermediateSteps,
            plotRhoA0Induced,
            saveIntermediatePlots,
            saveFinalPlots,
            saveIntermediatePlots
            )


