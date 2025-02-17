import numpy as np
import re
from pathlib import Path
from math import ceil

from plotScripts.readFiles import getDirectoryMA

floatRe = re.compile("\d+_\d+")

def floatToStr(value, sigDigs=3):
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

def strToFloat(value:str):
    return float(
        value.replace("_", ".")
        )

def openDataArray(
        dataArray: list,
        readData: bool=True,
        ):
    """
    This will be used to read all the Path objects in the readFiles dictionary.
    Parameters:
        dataArray, an array consisting of Path filenames.
        readData, if False, return the original array without doing anything
    Returns:
        dataArray, the read data array
    """
    if readData:
        dataArray = [
                np.genfromtxt(filename, delimiter=",")
                for filename in dataArray
                ]
    return np.array(dataArray)

def getLambdaValue(filename):
    # Each filename has only the mass and lambdaValue parameter
    return strToFloat(floatRe.findall(str(filename.name))[2])

def readSolutionsFromFile(
        m,
        a,
        lambdaValue,
        directory="",
        sigDigs=3,
        bcs="dirichlet",
        desiredQuantity:[str]=None,
        ):
    """
    Reads the converged solutions for a certain m, a, lambdaValue.
    """

    mString = floatToStr(m, sigDigs)
    lambdaString = floatToStr(lambdaValue, sigDigs)
    aString = floatToStr(a, sigDigs)

    fileId = f"mass_{mString}_a_{aString}_lambda_{lambdaString}.txt"

    if directory != "": 
        directory = "/" + directory

    # Initialize the solution family dict
    solutionFamily = {}

    for key in desiredQuantity:
        solutionFamily[key]  = np.genfromtxt(
            f"savedSolutions{directory}/{bcs}/{key}/{fileId}",
            dtype=float,
            delimiter="\n",
        )

    return solutionFamily

def downsizeUnreadSolutionFamily(
        posixDict,
        maxLambdaDensity,
        ):

    if maxLambdaDensity is None:
        return posixDict
    lambdaValueArray = posixDict["lambdaValue"]
    lambdaDensity =   len(lambdaValueArray) / (lambdaValueArray[-1] - lambdaValueArray[0])

    if lambdaDensity < maxLambdaDensity:
        # Do nothing if lambdaDensity is not too big,
        # or maxLambdaDensity was None
        return posixDict
    
    # Now lambdaDensity > maxLambdaDensity
    densityFactor = lambdaDensity / maxLambdaDensity
    # I want the len() of the arrays to be len(posixDict) / densityFactor
    newDict = {}
    for key, valueArray in posixDict.items():
        newDict[key] = valueArray[::ceil(densityFactor)]
    
    return newDict


def getPosixForQuantities(
        m,
        a,
        directory,
        sigDigs=3,
        maxLambdaDensity=None,
        ):


    dictKeys = ['eigenvalues', 'eigenstates',  'rho', 'A0Induced']
    # Not to fill the dictionary with crap that we don't need
    
    m = floatToStr(m, sigDigs=3)
    a = floatToStr(a, sigDigs=3)

    fileRegex = f"mass_{m}_a_{a}_lambda*"

    solutionFamilyDict = {}
    lambdaValueArray = []

    for i, quantity in enumerate(dictKeys):
        dataElement = sorted(
                list(
                    Path(f'{directory}/{quantity}').glob(fileRegex)
                    ),
                key = getLambdaValue
                )
        
        solutionFamilyDict[quantity] = dataElement

        # Enough with doing it once
        if i==0:
            lambdaValueArray = [
            getLambdaValue(filename) for filename in dataElement 
            ]

    solutionFamilyDict["lambdaValue"] = lambdaValueArray

    solutionFamilyDict = downsizeUnreadSolutionFamily(
            solutionFamilyDict, 
            maxLambdaDensity,
            )

    return solutionFamilyDict

def openPosixDict(
        posixDict,
        desiredQuantities: [str]=None,
        ):
    """
    Given a dictionary with families of solutions that store Path of solutions,
    open the desired quantities.
    Parameters:
        posixDict: A dictionary with the keys ['eigenvalues', 'eigenstates', 'eigenstateGradients', 'rho', 'A0Induced'],
                and items with the path of each of the files.
    """
    
    dictKeys = ['eigenvalues', 'eigenstates',  'rho', 'A0Induced']
    # Not to fill the dictionary with crap that we don't need
    finalDesiredQuantities = [] # To be able to ignore queries that were miswriten
    if desiredQuantities is None:
        finalDesiredQuantities = dictKeys
    else: 
        if isinstance(desiredQuantities, str):
            # Means desired quantities was a single str input
            desiredQuantities = [desiredQuantities]
        for quantity in desiredQuantities:
            if quantity not in dictKeys:
                print(f"Solicited quantity {quantity} is not accesible. I'm skipping it.")
            else:
                finalDesiredQuantities.append(quantity)

    solutionFamilyArray = {}

    for quantity in finalDesiredQuantities:
        dataElement = openDataArray(posixDict[quantity])
        solutionFamilyArray[quantity] = dataElement

    try:
        solutionFamilyArray["lambdaValue"] = posixDict["lambdaValue"]
    except KeyError:
        pass

    return solutionFamilyArray

def openSolutionFamilyArray(filterRegex) -> dict:
    """
    Selects and returns the solution family array, which is actually a dictionary, whose entries are arrays corresponding to each converged solution for each lambdaValue
    """

    directory, m, a = getDirectoryMA( filterRegex=filterRegex)
    posixDict = getPosixForQuantities(m, a, directory=directory)
    solutionFamilyArray = openPosixDict(posixDict)

    return directory, solutionFamilyArray

if __name__ == "__main__":
    m = 0
    a = 1
    lambdaValue = 7
    directory = "ambjorn"

    posixDict = getPosixForQuantities(m, a, directory=directory)
    
    solutionFamilyArray = openPosixDict(posixDict, desiredQuantities=["eigenvalues"])
    print(solutionFamilyArray)
