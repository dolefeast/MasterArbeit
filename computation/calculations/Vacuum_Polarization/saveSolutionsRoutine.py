import numpy as np
import os

def saveSolutions(
        self,
        sigDigs=3,
        ):
    """
    Saves the calculated quantities to {directory}/{boundaryConditions}/{quantity}/mass_{mass}A_{a}LambdaValue_{lambdaValue}.csv
    Parameters:
        solutionFamily: A dictionary with keys() = ["eigenvalues", "eigenstates", "eigenstateGradients"]
        directory: In case a further directory should be considered, e.g. if Ambjorn technique is used
    Returns None
    """

    solutionFamily = {
            "eigenvalues":self.eigenvalues,
            # "eigenstates":self.eigenstates,
            "A0Induced":self.A0Induced(self.z),
            "rho":self.rho,
            }

    if self.saveEigenstates:
        solutionFamily["eigenstates"] = self.eigenstates


    lambdaString = self.floatToStr(self.lambdaValue, sigDigs=sigDigs)
    aString = self.floatToStr(self.a, sigDigs=sigDigs)
    mString = self.floatToStr(self.m, sigDigs=sigDigs)

    fileId = f"mass_{mString}_a_{aString}_lambda_{lambdaString}.txt"

    if self.directory != "":
        directory = "/" + self.directory

    
    rootDirectory = f"savedSolutions{directory}/{self.bcs}"
    print(f"Saving results under {rootDirectory}/.../{fileId}...")
    
    try:
        for key, value in solutionFamily.items():
            np.savetxt(
                    f"{rootDirectory}/{key}/{fileId}",
                value,
                delimiter=",",
            )
    except FileNotFoundError as e:
        print(e)
        create = "y"
        if create == "y":
            for key in solutionFamily.keys():
                os.makedirs(rootDirectory + "/" + key)
            self.saveSolutions()
        elif create == "n": 
            rename = input(f"If {directory} was a typo, enter the correct name...")
            if rename != "":
                self.saveSolutions()
