import numpy as np
from scipy.interpolate import CubicSpline

def calculateRho(self, ):
    """
    Calculates the vacuum polarization AFTER normalizing the eigenstates, but before filtering
    """

    # The total charge density. 
    # Initiate as np.array to add them
    rho = np.zeros_like(self.z) 

    for n, (eigenvalue, eigenstate) in enumerate(
            zip(
                self.eigenvalues,
                self.eigenstates,
                )
            ):
        # The charge density associated to the nth mode
        # eigenstate is np array. A0 is a callable. 
        # A0(z) is of the same shape as eigenstate
        rhoN = (eigenvalue - self.A0(self.z)) * abs(eigenstate) ** 2
        if self.subtractMasslessPertVacuumPolarization and n>0:
            rhoN -=  self.perturbativeVacuumPolarizationNModeMassless(n, self.z) 
        rho += rhoN

    if self.subtractMasslessPertVacuumPolarization:
        rho +=  self.perturbativeVacuumPolarizationMasslessInf(self.z)

    return rho

def calculateDynamicRelaxParameter(self):
    """
    Given a previous A0 function of the system, calculate relaxParameter 
    as the value for which the derivative of the update rule f(A0) is 0.

    For this calculate "The next A0" for two different given potentials.
    In this case, they stem from the same potential, but are both slightly perturbed.
    """ 

    A1Induced = lambda z: 1.5 * self.A0History[-1][-1](z) 
    
    self.A0 = lambda z: -self.lambdaValue * (z-1/2) + A1Induced(z)
    self.A0Induced = A1Induced
    self.coreIteration()
    A1InducedNew = self.A0Induced

    A2Induced = lambda z: 0.9 * self.A0History[-1][-1](z) 
    self.A0 = lambda z: -self.lambdaValue * (z-1/2) + A2Induced(z)
    self.A0Induced = A2Induced
    self.coreIteration()
    A2InducedNew = self.A0Induced

    c = lambda z: (A1InducedNew(z) - A2InducedNew(z)) / (A1Induced(z) - A2Induced(z) - A1InducedNew(z) + A2InducedNew(z))
    print(c(1))
    return c

def calculateRelaxParameter(self):
    """
    In the update law A_{k+1} = (1-c) AK + c ( - λ (z - 1/2) - ∫ ∫ ρ ), 
    calculate c so that the correction term c( - λ ... ) is exactly self.correctionParameter
    """

    if self.relaxParameter is None:
        c = (
                self.correctionParameter * (- self.lambdaValue / 2 + self.a * self.A0Induced(1) )  
                / ( (1-self.correctionParameter) * self.constantLambdaA0List[-1](1) 
                    + self.correctionParameter * ( - self.lambdaValue / 2 + self.a * self.A0Induced(1))
                    )
                )
    elif self.dynamicRelaxParameter:
        c = self.calculateDynamicRelaxParameter()
    else: 
        c = self.relaxParameter
    
    return c

def calculateA0Induced(self, rho):
    """
    Calculates the induced A0 as the solution to the Poisson equation.
    It is just integrating the vacuum polarization twice
    """

    rhoInterpolated = CubicSpline(self.z, rho)

    A0InducedShifted = rhoInterpolated.antiderivative(2)
    offset = A0InducedShifted(1/2)
    A0Induced = lambda z: -(A0InducedShifted(z) - offset)

    if A0Induced(1) > 10:
        print(f"n={self.n}, A0Induced(1)={self.A0Induced(1)}")

    return A0Induced

def initializeA0Lists(self, extrapolate=False):
    
    # Need at least two converged potentials for this scheme to work
    extrapolate = extrapolate and len(self.lambdaArray) > 2
    # Calculate the new A0 by point-wise linear extrapolation from the last two A0
    if extrapolate:
        # The last lambda value, and the one before that
        lambda1 = self.lambdaArray[-1] 
        lambda2 = self.lambdaArray[-2]

        # The last full potential A0, and the one before that
        A1 = self.A0History[-1][-1]
        A2 = self.A0History[-2][-1]

        # The last induced potential A0, and the one before that
        A0Induced1 = self.A0InducedHistory[-1][-1]
        A0Induced2 = self.A0InducedHistory[-2][-1]

        A0 = lambda z: (A1(z)- A2(z))/(lambda1 - lambda2) * (self.lambdaValue - lambda1) + A1(z)
        A0Induced = lambda z: (A0Induced1(z)- A0Induced2(z))/(lambda1 - lambda2) * (self.lambdaValue - lambda1) + A0Induced1(z)
    else:
        A0 = self.A0History[-1][-1] 
        A0Induced =  self.A0InducedHistory[-1][-1] 

    self.constantLambdaA0List = [ A0 ]
    self.constantLambdaA0InducedList = [ A0Induced ]

