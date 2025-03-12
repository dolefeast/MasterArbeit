# Some util functions
def sign(self, x):
    # Returns sign of x
    return 2 * (x>0) - 1 

def floatToStr(self, value, sigDigs=3):
    if isinstance(value, str):
        return value
    return str(
        round(
            float(
                value
            ),
            sigDigs
        )
    ).replace(".", "_")

def strToFloat(self, value:str):
    return float(
        value.replace("_", ".")
        )

def rootMeanSquare(self, y):
    return np.sqrt(np.mean(y**2))

def setConfigFromDict(self, configDict:dict, calculateEigenstates: bool, modifyParameter:float=1):
    from scipy.interpolate import CubicSpline
    from numpy import linspace

    self.eigenvalues = configDict["eigenvalues"]
    self.maxN = len(self.eigenvalues)//2
    if self.bcs == "neumann":
        self.maxN -= 1 # Because neumann includes a zeroth case

    self.nPoints = len(configDict["A0Induced"])
    self.z = linspace(0, 1, self.nPoints)

    # These quantities here are really only for the user's convenience
    self.A0Induced = CubicSpline(self.z, modifyParameter * configDict["A0Induced"] )
    self.A0 = lambda z: - self.lambdaValue * (z - 1/2) + self.e * self.A0Induced(z)
    self.rho =  configDict["rho"]
    self.lambdaValue = configDict["lambdaValue"]

    # This is what actually gets used. Avoids recursion problems
    self.A0InducedHistory = [[ self.A0Induced ]]
    self.A0History = [[ self.A0 ]]

    if calculateEigenstates:
        self.calculateEigenstatesParallel()
