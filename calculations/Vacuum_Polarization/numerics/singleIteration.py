def singleIteration(self):

    print("Calculating eigenstates")
    self.calculateEigenstatesParallel()

    print("Calculating rho")
    self.calculateRho()
    
    if self.smoothing:
        print("Smoothing")
        self.convolveRho()

    self.calculateA0Induced()
