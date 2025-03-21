import numpy as np
from pathlib import Path, PurePath
import re

floatRe = re.compile("\d+_\d+")

def getDirectoryMA(readMode='str', filterRegex=""):
    dirs = Path("data")

    if filterRegex != "":
        filterRegex = "*" + filterRegex + "*"
    else:
        filterRegex = "*"

    [print(f"[{i+1}]", dirName.name) for i, dirName in enumerate(sorted(list(dirs.glob(filterRegex))))]

    while True:
        try:
            indexSTR = input("Choose the desired directory... (empty input to stop)")
            if indexSTR == "":
                return 
            index = int(indexSTR)
            directory = sorted(list(dirs.glob(filterRegex)))[index-1]
            break
        except ValueError:
            print("Input must be a natural number, try again")
        except IndexError:
            print(f"Input must be a number from 1 to {len(list(dirs.glob('*')))}")

    bcsList = list(directory.glob("*"))

    if len(bcsList) == 1:
        bcs = str(bcsList[0].name)
        print(f"The only boundary conditions found are {bcs}")
    else:
        print("\nAvailable boundary conditions are")
        [print(f"[{i+1}]", dirName.name) for i, dirName in enumerate(bcsList)]
        while True:
            try:
                indexBcs = int(input("Choose the desired boundary conditions... "))
                bcs = str(bcsList[indexBcs-1].name)
                break
            except ValueError:
                print("Input must be a natural number, try again")
            except IndexError:
                print(f"Input must be a number from 1 to {len(bcsList)}")
        print(f"\nChosen boundary conditions: {bcs}")


    # Since the namefiles are repeated in the different quantity names I only need to do this once 
    # to get all different m, a values.
    for name in ["rho"]: 
        quantityDir = Path(PurePath(directory) / bcs / name)

        mASet = set()
        for filename in quantityDir.glob("**/*txt"):
            try:
                if readMode=='str':
                    matches = re.findall("\d+_\d+", str(filename.name))
                    m, a, lambdaValue = matches
                elif readMode=='float':
                    m, a, lambdaValue = map(strToFloat, re.findall("\d+_\d+", str(filename.name)))
                else:
                    raise ValueError(f"readMode must be either 'str', 'float' but it was {readMode}")
            except ValueError as e:
                print(e)
                continue
            mASet.add((m,a))

    if len(mASet) > 1:
        print("The distinct m, a values are:")
        for i, mA in enumerate(mASet):
            print(f"[{i+1}]: m = {list(mA)[0]}, a = {list(mA)[1]}") 
        while True:
            try:
                index = int(input("Choose the desired m, a value ... "))
                mAChoice = list(mASet)[index-1]
                break
            except ValueError:
                print("Input must be a natural number, try again")
            except IndexError:
                print(f"Input must be a number from 1 to {len(list(dirs.glob('*')))}")
    else:
        print(f"There is only one distinct case (m, a) = {list(mASet)[0]}")
        mAChoice = list(mASet)[0]

    return directory/bcs, mAChoice[0], mAChoice[1]
def getIndexForValueFromBelow(array, value):
    """
    Given an array, return the biggest index i for which array[i] <= value
    """
    for i, element in enumerate(array):
        if element > value:
            return i-1
    return -1

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

def getLambdaValue(self,filename):
    # Each filename has only the mass and lambdaValue parameter
    return self.strToFloat(floatRe.findall(str(filename.name))[2])

def readSolutionsFromFile(
        self,
        m,
        a,
        lambdaValue,
        directory="",
        bcs="dirichlet",
        desiredQuantity:[str]=None,
        ):
    """
    Reads the converged solutions for a certain m, a, lambdaValue.
    """
    sigDigs = self.sigDigs

    mString = self.floatToStr(m, sigDigs)
    lambdaString = self.floatToStr(lambdaValue, sigDigs)
    aString = self.floatToStr(a, sigDigs)

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
        self,
        m,
        a,
        directory,
        sigDigs=3,
        maxLambdaDensity=None,
        ):


    dictKeys = ['eigenvalues', 'eigenstates',  'rho', 'A0Induced']
    # Not to fill the dictionary with crap that we don't need
    
    m = self.floatToStr(m)
    a = self.floatToStr(a)

    fileRegex = f"mass_{m}_a_{a}_lambda*"

    solutionFamilyDict = {}
    lambdaValueArray = []

    for i, quantity in enumerate(dictKeys):
        dataElement = sorted(
                list(
                    Path(f'{directory}/{quantity}').glob(fileRegex)
                    ),
                key = self.getLambdaValue
                )
        
        solutionFamilyDict[quantity] = dataElement

        # Enough with doing it once
        if i==0:
            lambdaValueArray = [
            self.getLambdaValue(filename) for filename in dataElement 
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
    
    dictKeys = ['eigenvalues', 'rho', 'A0Induced', "eigenstates"]
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

def dictIndex(dictionary, index):
    for key, value in dictionary.items():
        try:
            dictionary[key] = value[index]
        except IndexError:  # The desired index was too big
            try:
                dictionary[key] = value[-1]
            except IndexError: # The array is empty
                pass

    return dictionary

def openSolutionFamilyArray(self, filterRegex, index=None) -> dict:
    """
    Selects and returns the solution family array, which is actually a dictionary, whose entries are arrays corresponding to each converged solution for each lambdaValue
    """

    directory, m, a = getDirectoryMA( filterRegex=filterRegex)
    posixDict = self.getPosixForQuantities(m, a, directory=directory)
    solutionFamilyArray = openPosixDict(posixDict)

    if index is not None:
        return directory, self.strToFloat(m) , self.strToFloat(a) , dictIndex(solutionFamilyArray, index)

    return directory, self.strToFloat(m) , self.strToFloat(a) , solutionFamilyArray


def setConfigFromDict(self, configDict:dict=None, calculateEigenstates=False, modifyParameter:float=1, index=-1):
    from scipy.interpolate import CubicSpline
    from numpy import linspace

    if configDict is None:
        directory, self.m, self.a, configDict = self.openSolutionFamilyArray(filterRegex="", index=index)

    self.eigenvalues = configDict["eigenvalues"]
    self.maxN = len(self.eigenvalues)
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
