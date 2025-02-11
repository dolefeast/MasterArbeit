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
            lambdaMin=10,
            lambdaMax=25,
            lambdaStep=2,
            relaxParameter=0.2,
            maxNIterations=200,
            walkbackFactor=0.4,
            lambdaStepMin=1e-6,
            nIterations=None,
            iterateTol=1e-3,
            smoothing=False,
            ambjorn=False,
            plotForEachLambda = False, 
            saveData = False, 
            savePlots = False, 
            showPlots = False, 
            read = False, 
            directory = "",
            bcs="dirichlet",
            extrapolateA0=False,
            dynamicRelaxParameter=False,
            subtractMasslessPertVacuumPolarization=False,
            ):

        self.e = e
        self.a = a
        self.maxN = maxN
        self.m = m
        self.eigenvalues = [ 
            self.perturbativeEigenvalue(n) for n in range(-maxN, maxN+1) if n!=0
            ]

        if bcs == "neumann":
            self.eigenvalues = self.eigenvalues[:self.maxN] + [-self.m, self.m] + self.eigenvalues[self.maxN:]


        self.colorCycle = cycle(
    "#C52E19FF, #AC9765FF, #54D8B1FF, #B67C3BFF, #175149FF, #AF4E24FF".split(", ")
        )
        self.color = next(self.colorCycle)

         
        self.plot = savePlots or showPlots
        self.plotForEachLambda = plotForEachLambda
        self.saveData = saveData
        self.savePlots = savePlots
        self.showPlots = showPlots
        self.read = read

        if self.plot:
            import matplotlib.pyplot as plt
            self.plt = plt
            self.plottingSetup()

        self.relaxParameter = relaxParameter # the number c in A_{k+1} = (1-c) AK + c ( -𝜆 (z-1/2) - ∫∫𝜌)
        self.extrapolateA0 = extrapolateA0
        self.dynamicRelaxParameter = dynamicRelaxParameter # Whether c gets calculated at each lambda value to try to make things go faster

        if maxN<5:
            self.nPoints = 200
            self.ambjorn = True # This can't be changed
            self.smoothing = False
        else:
            self.nPoints = 8 * ( maxN + 1 )
            self.ambjorn = ambjorn      # This can be changed depending on the results we are looking for
            self.smoothing = smoothing  # This can be changed depending on the results we are looking for
        

        # Part of the vacuum polarization filtering routine
        self.subtractMasslessPertVacuumPolarization = self.smoothing and subtractMasslessPertVacuumPolarization

        # The mesh points
        self.z = np.linspace(0, 1, self.nPoints)
        
        self.lambdaMin = lambdaMin
        self.lambdaMax = lambdaMax
        self.lambdaStep = lambdaStep
        self.lambdaValue = lambdaMin
        self.lambdaArray = []

        self.rho = [0] * self.nPoints
        self.A0Induced = lambda z: 0
        self.A0InducedHistory = [[self.A0Induced]]
        self.A0 = lambda z: - self.lambdaValue * (z - 1/2) 
        self.A0History = [[self.A0]]

        self.nIterations = nIterations # In case I want it to just iterate a small finite amount of times
                    # if nIterations = None, "back-react" until convergence
        self.iterateTol = iterateTol # Tolerance for the convergence of the backreaction procedure
        self.maxNIterations = maxNIterations # Avoid the calculations from getting in a loop
        self.walkbackFactor = walkbackFactor # Avoid the calculations from getting in a loop
        self.lambdaStepMin = lambdaStepMin # Avoid the calculations from getting in a loop

        self.directory = directory
        self.bcs = bcs

        if savePlots:
            directoryExists = Path("figures/"+directory+"/"+bcs).is_dir()

            if directoryExists:
                overwrite = "y"
                if overwrite=="y":
                    shutil.rmtree("figures/"+directory+"/"+bcs)
                elif overwrite=="n":
                    print("Then change the code")
                    exit()

            Path("figures/"+directory+"/"+bcs).mkdir(parents=True, exist_ok=True)
    from Vacuum_Polarization.filterRoutines import (
            extendSignal, 
            removeNeighbourhood,
            removeAndInterpolate,
            returnTo01,
            filterRho,
            extendAndFilter 
            )

    from Vacuum_Polarization.smallRoutines import (
            sign,
            floatToStr,
            strToFloat,
            rootMeanSquare,
            setConfigFromDict
            )

    from Vacuum_Polarization.plottingSetup import (
            plottingSetup,
            plotRhoA0Induced,
            plotIntermediateStepsSetup,
            plotIntermediateSteps,
            saveAllPlots
            )
    from Vacuum_Polarization.bisectionRoutines import (
            findRootBisection,
            bisectionMethodUpperBound,
            bisectionMethodLowerBound,
            )
    from Vacuum_Polarization.physics import (
            KleinGordonEquation
            )
    from Vacuum_Polarization.perturbativeSolutions import(
            perturbativeEigenvalue,
            perturbativePhi,
            freePhi,
            perturbativeVacuumPolarizationMasslessInf,
            perturbativeVacuumPolarizationNModeMassless,
            )

    from Vacuum_Polarization.calculateEigenstatesRoutines import (
            calculateEigenstates,
            normalizeEigenstates,
            )

    from Vacuum_Polarization.saveSolutionsRoutine import saveSolutions

    from Vacuum_Polarization.inducedFieldsRoutines import (
            calculateRho,
            calculateRelaxParameter,
            calculateA0Induced,
            initializeA0Lists,
            calculateDynamicRelaxParameter
            )
    from Vacuum_Polarization.mainScript import (
            coreIteration,
            constantLambdaIterations,
            walkback,
            fullScript,
            )
