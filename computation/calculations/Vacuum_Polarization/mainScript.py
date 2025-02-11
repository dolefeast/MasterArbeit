from scipy.interpolate import CubicSpline
import numpy as np

def coreIteration(self):
    """
    Calculates and normalizes the eigenstates for the current A0.
    Calculates the vacuum polarization corresponding to this solution family, 
    filters it if necessary, and finally calculates the corresponding induced potential to this A0Induced.
    """
    if self.parallelization:
        self.calculateEigenstatesParallel()
    else:
        self.calculateEigenstates()
        self.normalizeEigenstates()

    rhoNew = self.calculateRho()

    if not self.ambjorn:
        rhoNew += self.e**2/np.pi * self.A0(self.z)

    if self.smoothing:
        # rho = filterRho(rho)
        _, rhoNew = self.extendAndFilter(rhoNew)

    if False and self.n>70 and not self.n%25:
        print("Averaging..")
        self.rho = (rhoNew + self.rho)/2 
    else:
        self.rho = rhoNew

    self.A0Induced = self.calculateA0Induced(self.rho)

def constantLambdaIterations(self):
    """
    Calculates until convergence or until a certain number of iterations, the backreaction iteration.
    During this calculation it might save different plots for each of the iterations.
    """
    # Checks the convergence of rho as the iterations go
    rhoError=1

    if self.savePlots and self.plotForEachLambda:
        self.plotIntermediateStepsSetup()

    self.n=0

    self.initializeA0Lists(extrapolate=self.extrapolateA0)

    while self.n < self.maxNIterations:
        c = self.calculateRelaxParameter()

        self.A0Induced = self.constantLambdaA0InducedList[-1]
        if callable(c):
            self.A0 = lambda z:  (
                    (1-c(z)) * self.constantLambdaA0List[-1](z) + c(z) * (- self.lambdaValue * (z-1/2) + self.a * self.A0Induced(z)) 
                    )
        else:
            self.A0 = lambda z:  ( (1-c) * self.constantLambdaA0List[-1](z) + c * (- self.lambdaValue * (z-1/2) + self.a * self.A0Induced(z)) ) 

        # If we are only looking for a fixed amount of iterations:
        if isinstance(self.nIterations, int):
            if self.n >= self.nIterations: break
        # Or if we want to look for convergence (i.e. nIterations -> infty )
        elif self.nIterations is None:
            if rhoError < self.iterateTol: 
                if convergedCount == 2:
                    break
                else:
                    convergedCount += 1
            else:
                convergedCount = 0
        else: 
            raise Exception("nIterations must be either int or None")

        # The main calculation. Calculate KG solutions, the associated rho and the associated A0Induced
        self.coreIteration()
        
        # If we want to plot things
        self.plotIntermediateSteps() 

        if self.nIterations is None:
            try:
                # Tries to calculate the error if maxPrevRho is defined
                rhoError = abs(max(self.rho) - maxPrevRho) / abs(maxPrevRho)
            except NameError:
                # If it wasn't, define it and calculate again
                pass
            maxPrevRho = max(self.rho)

        self.A0 = CubicSpline(self.z, self.A0(self.z))

        self.constantLambdaA0List.append(self.A0)
        self.constantLambdaA0InducedList.append(self.A0Induced)
        self.n += 1

        print(self.n, end = "\r")
    else:
        if self.savePlots:
            self.fig4.savefig(
                    "figures/"+self.directory+"/"+self.bcs+"/intermediateStepsLambdaValue_"+self.floatToStr(self.lambdaValue)+".png"
                    )
        raise RecursionError("maxNiterations reached")

    if self.savePlots: 
        self.fig4.savefig(
                "figures/"+self.directory+"/"+self.bcs+"/intermediateStepsLambdaValue_"+self.floatToStr(self.lambdaValue)+".png")

def walkback(self):
    """
    If a failure in the calculation was found, "walk back" i.e. return to the last lambdaValue that lead to converged solutions, and reduce 1. the lambdaStep 2. the relax parameter in the update law
    """

    self.lambdaValue -= self.lambdaStep 
    self.lambdaStep *= self.walkbackFactor
    self.lambdaValue += self.lambdaStep
    if self.a!=0:
        self.relaxParameter *= self.walkbackFactor
        print(f"New relaxParameter={self.relaxParameter}")
def fullScript(self):

    while self.lambdaValue < self.lambdaMax:
        if self.lambdaStep < self.lambdaStepMin:
            e = Exception("lambdaStep got too small")
            break

        print(15*"#")
        print('λ =', self.lambdaValue)
        self.color = next(self.colorCycle) # To distinguish between different lambda values

        try:
            self.constantLambdaIterations() 

            # If no exception was raised, a self-consistent solution was found and it converged
            self.lambdaArray.append(self.lambdaValue)

            if self.plot:
                self.eigenvaluesArray.append(self.eigenvalues)
                self.plotRhoA0Induced(self.rho, self.A0Induced(self.z), '', color=self.color, alpha=1)

            print(f'Converged in n={self.n} iterations. A0Induced(1) =', self.A0Induced(1))
            if self.saveData:
                self.saveSolutions()

        except RecursionError as exception:
            self.walkback()

            print(exception) # To know if it didn't converge or if no solution was found
            print(f'Error found at n={self.n}, A0Induced(1) =', self.A0Induced(1))
            continue
        except RuntimeError as exception:
            self.walkback()

            print(exception) # To know if it didn't converge or if no solution was found
            print(f'Error found at n={self.n}, A0Induced(1) =', self.A0Induced(1))
            continue
        except ValueError as exception:
            e = exception
            break
        except Exception as exception:
            e = exception
            break
        except KeyboardInterrupt as exception:
            e = exception
            break

        self.A0History.append(self.constantLambdaA0List)
        self.A0InducedHistory.append(self.constantLambdaA0InducedList)


        self.lambdaValue += self.lambdaStep
    
    if self.savePlots:
        self.saveAllPlots()
    if self.showPlots:
        self.axEigenvalues.plot(self.lambdaArray, self.eigenvaluesArray, 'b')
        self.plt.show()
    try:
        raise e 
    except UnboundLocalError:
        print("Routine reached lambdaMax succesfully")
