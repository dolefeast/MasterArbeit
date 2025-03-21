def initializeA0Lists(self):
    # Need at least two converged potentials for this scheme to work
    extrapolate = self.extrapolate and len(self.lambdaArray) > 2
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

def convergence(self):
    from scipy.interpolate import CubicSpline
    """
    Calcualtes new Klein-Gordon solutions until convergence
    """

    # Initializes the lists self.constantLambdaA0List and self.constantLambdaA0InducedList
    # Also tries to predict the next A0 via extrapolation
    self.initializeA0Lists()  

    self.n = 0

    while self.n <= self.maxNIterations:

        self.A0Induced = self.constantLambdaA0InducedList[-1]
        c = self.relaxParameter # For the sake of brevity
        
        self.A0 = lambda z:  ( 
                (1-c) * self.constantLambdaA0List[-1](z) 
                + c * (- self.lambdaValue * (z-1/2) + self.a * self.A0Induced(z)) 
                )

        if self.nIterations is None:
            if self.bcs == "dirichlet" and self.m == 0 and self.lambdaValue < 10 and self.n==1:
                # Very special cases that don't require iterations
                break
        else:
            if self.n >= self.nIterations: break

        try:
            self.singleIteration()
        except self.NoRootFoundError as e:
            print(e)
            return False
        except self.MaxIterError as e:
            print(e)
            return False


        
        self.plotIntermediateSteps() # Plots the intermediate values 
                                                           # Allows to see convergence, etc.
        if self.nIterations is None:
            try:
                # Tries to calculate the error if maxPrevRho is defined
                rhoError = abs(max(self.rho) - maxPrevRho) / abs(maxPrevRho)
                if rhoError < self.convergenceTol:
                    break # Solution converged
            except NameError:
                # If it wasn't, define it and calculate again
                pass
            maxPrevRho = max(self.rho)

        self.A0 = CubicSpline(self.z, self.A0(self.z)) # Fixes it and avoids recursion errors.

        self.constantLambdaA0List.append(self.A0)
        self.constantLambdaA0InducedList.append(self.A0Induced)
        self.n += 1

    else: # Converged too slowly
        self.saveIntermediatePlots()  
        return False

    self.saveIntermediatePlots()  
    return True

