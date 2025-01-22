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

def setConfigFromDict(self, configDict:dict):
    from scipy.interpolate import CubicSpline
    from numpy import linspace

    self.eigenvalues = configDict["eigenvalues"]
    self.maxN = len(self.eigenvalues)//2
    self.eigenstates = configDict["eigenstates"]
    self.nPoints = len(configDict["A0Induced"])
    self.z = linspace(0, 1, self.nPoints)

    self.A0Induced = CubicSpline(self.z, configDict["A0Induced"])
    self.rho =  configDict["rho"]
    self.lambdaValue = configDict["lambdaValue"]
