from pathlib import Path, PurePath
import re
import matplotlib.pyplot as plt
<<<<<<< HEAD

import matplotlib_format
=======
import scienceplots

plt.style.use(["science", "high-contrast"])
plt.rcParams["figure.figsize"] = (3.5, 2.6)
plt.rcParams["font.size"] = "5.4"
plt.rcParams["axes.labelsize"] = "13"
plt.rcParams["xtick.labelsize"] = "13"
plt.rcParams["ytick.labelsize"] = "13"
plt.rcParams["lines.linewidth"] = "0.9"
>>>>>>> 56f278ba38b8e6dca00bf7b5f466caed29774e7e

thingsRe = re.compile(r"([a-zA-Z0]*)_lambda_([\d.]*)_([a-z]+)")

dirs = Path("convergenceLog")

figOmega, axOmega = plt.subplots()
figA0, axA0 = plt.subplots()


figDict = {
            "eigenvalues":figOmega,
            "A0Induced":figA0
          }

axDict = {
            "eigenvalues":axOmega,
            "A0Induced":axA0
          }

axTitle = {
    "eigenvalues": r"$\omega^\kappa_\lambda$",
    "A0Induced": r"$A_{0, \lambda}^{\text{br},\kappa}(1)$",
}

for filename in list(dirs.glob("*"))[:]:
    quantity, lambdaValue, bcs = re.findall(thingsRe, str(filename))[0] 
    lambdaValue = round(float(lambdaValue), 3)
    
    kappa = 1
    ax = axDict[quantity]
    fig = figDict[quantity]
    fig.suptitle(f"$\lambda = {lambdaValue}$")
    ax.cla()
    ax.set_xlabel("$\kappa$")
    ax.set_ylabel(axTitle[quantity])

    with open(filename, "r") as myFile:
        for line in myFile:
            if quantity == "A0Induced":
                inv = -1
            else:
                inv = 1
            ax.scatter(kappa, inv*float(line.split(", ")[0]), color="#0E103D")
            kappa += 1
    
    fig.tight_layout()
    fig.savefig(f'convergenceLogFigs/{filename.stem}.pdf')

