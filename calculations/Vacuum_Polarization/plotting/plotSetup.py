import matplotlib.pyplot as plt
def plotSetup(self):
    self.figRho = plt.figure(r"Vacuum polarization")
    self.figA0Induced = plt.figure(r"A_0 induced")
    self.figEigenvalues = plt.figure("Mode energy evolution")

    if self.saveIntermediatePlots:
        self.figIntermediateA0 = plt.figure()
        self.figIntermediateEigenvalue = plt.figure()

    self.axRho = self.figRho.subplots(1)
    self.axA0Induced = self.figA0Induced.subplots(1)
    self.axEigenvalues = self.figEigenvalues.subplots(1)

    self.axRho.set_ylabel(r"$\rho$")
    self.axRho.set_xlabel(r"z")

    self.axA0Induced.set_ylabel(r"$A_0^{br}$")
    self.axA0Induced.set_xlabel(r"z")

    self.axEigenvalues.set_ylabel(r"$\omega_n$")
    self.axEigenvalues.set_xlabel(r"$\lambda$")

    self.eigenvaluesArray = []

def plotIntermediateStepsSetup(self):
    self.figIntermediateA0.suptitle(r"$\lambda$={}".format(round(self.lambdaValue, 4)), fontsize=10)
    self.figIntermediateEigenvalues.suptitle(r"$\lambda$={}".format(round(self.lambdaValue, 4)), fontsize=10)

    self.axRelax.cla()
    self.axIntermediateOmegas.cla()

    self.axIntermediateOmegas.set_xlabel(r"Iteration")
    self.axRelax.set_xlabel(r"Iteration")
    self.axRelax.set_ylabel(r"$A0(1) +$ {}".format(round(self.lambdaValue/2, 5)))

    # The mode we are most interested in when studying the different boundary conditions
    if self.bcs == "dirichlet":
        mode = "1"
    elif self.bcs == "neumann":
        mode = "0"

    self.axIntermediateOmegas.set_ylabel(r"$\omega_{}$".format(mode))
    self.axIntermediateOmegas.set_xlabel(r"Iteration".format(mode))

def plotRhoA0Induced(self, fmt, color, alpha,):
    self.axRho.plot(self.z, self.rho, fmt, color=self.color, alpha=alpha)
    self.axA0Induced.plot(self.z, self.A0Induced, fmt, color=self.color, alpha=alpha)

def plotIntermediateSteps(self):
    if not self.intermediatePlots:
        return

    self.plotRhoA0Induced(self.rho, self.A0Induced(self.z), '-', color=self.color, alpha=0.2)

    color = "#0A2463"
    self.axRelax.plot(self.n, self.A0(1) + self.lambdaValue/2, 'o', color=color)
    # self.axRelax.plot(self.n, self.A0(1) + round(self.lambdaValue/2, 2), 'o', color=self.color)
    self.axIntermediateOmegas.plot(self.n, self.eigenvalues[0], 'o', color=color)

def saveFinalPlots(self):
    self.axEigenvalues.plot(self.lambdaArray, self.eigenvaluesArray, 'b')
    self.figRho.savefig("figures/"+self.directory+"/"+self.bcs+"/vacuumPolarization.pdf")
    self.figA0Induced.savefig("figures/"+self.directory+"/"+self.bcs+"/A0Induced.pdf")
    self.figEigenvalues.savefig("figures/"+self.directory+"/"+self.bcs+"/eigenvalues.pdf")

def saveIntermediatePlots(self):
    if not self.savePlots:
        return 
    self.figIntermediateA0.savefig(
            f"figures/{self.directory}/{self.bcs}/A0InducedintermediateStepsLambdaValue_{self.floatToStr(self.lambdaValue)}.pdf")
    self.figIntermediateEigenvalues.savefig(
            f"figures/{self.directory}/{self.bcs}/omegaintermediateStepsLambdaValue_{self.floatToStr(self.lambdaValue)}.pdf")
            
