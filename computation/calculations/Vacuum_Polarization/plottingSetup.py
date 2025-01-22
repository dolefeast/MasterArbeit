import matplotlib.pyplot as plt

def plottingSetup(self):
    self.fig1 = self.plt.figure(r"Vacuum polarization")
    self.fig2 = self.plt.figure(r"A_0 induced")
    self.fig3 = self.plt.figure("Mode energy evolution")
    self.fig4 = self.plt.figure("Values at intermediate steps")

    self.axRho = self.fig1.subplots(1)
    self.axA0Induced = self.fig2.subplots(1)
    self.axEigenvalues = self.fig3.subplots(1)
    self.axRelax, self.axIntermediateOmegas = self.fig4.subplots(2)

    self.axRho.set_ylabel(r"$\rho$")
    self.axRho.set_xlabel(r"z")

    self.axA0Induced.set_ylabel(r"$A_0^{br}$")
    self.axA0Induced.set_xlabel(r"z")

    self.axEigenvalues.set_ylabel(r"$\omega_n$")
    self.axEigenvalues.set_xlabel(r"$\lambda$")

    self.eigenvaluesArray = []

def plotRhoA0Induced(self, rho, A0Induced, fmt, color, alpha,):
    self.axRho.plot(self.z, rho, fmt, color=self.color, alpha=alpha)
    self.axA0Induced.plot(self.z, A0Induced, fmt, color=self.color, alpha=alpha)

def plotIntermediateStepsSetup(self):
    self.fig4.suptitle(r"$\lambda$={}".format(round(self.lambdaValue, 4)))

    self.axRelax.cla()
    self.axIntermediateOmegas.cla()

    self.axIntermediateOmegas.set_xlabel(r"n")
    self.axRelax.set_ylabel(r"$A0(1) +$ {}".format(round(self.lambdaValue/2, 5)))
    self.axIntermediateOmegas.set_ylabel(r"$\omega$")

def plotIntermediateSteps(self):
    if self.plot:
        self.plotRhoA0Induced(self.rho, self.A0Induced(self.z), '-', color=self.color, alpha=0.2)

    if self.plotForEachLambda and self.savePlots:
        self.axRelax.plot(self.n, self.A0(1) + round(self.lambdaValue/2, 2), 'x', color=self.color)
        self.axIntermediateOmegas.plot(self.n, self.eigenvalues[self.maxN], 'x', color=self.color)


def saveAllPlots(self):
    self.axEigenvalues.plot(self.lambdaArray, self.eigenvaluesArray, 'b')
    self.fig1.savefig("figures/"+self.directory+"/"+self.bcs+"/vacuumPolarization.png")
    self.fig2.savefig("figures/"+self.directory+"/"+self.bcs+"/A0Induced.png")
    self.fig3.savefig("figures/"+self.directory+"/"+self.bcs+"/eigenvalues.png")
