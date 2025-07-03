import numpy as np
import os

def saveDataScript(self):

    solutionFamily = {
            "eigenvalues":self.eigenvalues,
            "eigenstates":self.eigenstates[0], # This is only for visualization purposes
            "A0Induced":self.A0Induced(self.z),
            "rho":self.rho,
            }

    lambdaString = self.floatToStr(self.lambdaValue)
    aString = self.floatToStr(self.a)
    mString = self.floatToStr(self.m)

    fileId = f"mass_{mString}_a_{aString}_lambda_{lambdaString}.txt"
    
    rootDirectory = f"data/{self.directory}/{self.bcs}"
    try:
        for key, value in solutionFamily.items():
            np.savetxt(
                    f"{rootDirectory}/{key}/{fileId}",
                value,
                delimiter=",",
                fmt="%.5f"
            )
        print(f"Saving results under {rootDirectory}/.../{fileId}...")
    except FileNotFoundError:
        for key in solutionFamily.keys():
            os.makedirs(f"{rootDirectory}/{key}")
        self.saveDataScript()
