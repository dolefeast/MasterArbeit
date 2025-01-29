from pathlib import Path, PurePath
import re

quantityNames = ['eigenvalues', 'eigenstates',  'rho', 'A0Induced']
def getDirectoryMA(readMode='str', filterRegex=""):
    dirs = Path("savedSolutions")

    if filterRegex != "":
        filterRegex = "*" + filterRegex + "*"
    else:
        filterRegex = "*"

    [print(f"[{i+1}]", dirName.name) for i, dirName in enumerate(dirs.glob(filterRegex))]

    while True:
        try:
            index = int(input("Choose the desired directory... "))
            directory = list(dirs.glob(filterRegex))[index-1]
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

if __name__ == "_Main__":
    bcs = "dirichlet"
    print(getDirectoryName(bcs))
