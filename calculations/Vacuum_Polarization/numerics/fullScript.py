def fullScript(self):
    while self.lambdaValue < self.lambdaMax:
        if self.lambdaStep < self.lambdaStepMin:
            print("Lambda step got too small")
            return False
        print("\n", 15*"#")
        print('λ =', self.lambdaValue)

        # To distinguish between different lambda values
        self.color = next(self.colorCycle) 
        
        
        try:
            converged = self.convergence()
        except KeyboardInterrupt:
            print("que pasooooooooooooo")
            break

        if not converged:
            print(f'Error found at n={self.n}, A0Induced(1) =', self.A0Induced(1))
            self.walkback()
            continue

        if self.saveData:
            self.saveDataScript()

        self.A0History.append(self.constantLambdaA0List)
        self.A0InducedHistory.append(self.constantLambdaA0InducedList)

        self.lambdaValue += self.lambdaStep

    if self.showPlots:
        self.axEigenvalues.plot(self.lambdaArray, self.eigenvaluesArray, 'b')
        plt.show()
    if self.savePlots:
        self.saveAllPlots()


def walkback(self):
    self.lambdaValue -= self.lambdaStep 
    self.lambdaStep *= self.walkbackFactor
    self.lambdaValue += self.lambdaStep

    if self.a!=0:
        self.relaxParameter *= self.walkbackFactor

