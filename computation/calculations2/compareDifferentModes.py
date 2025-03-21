import matplotlib.pyplot as plt
import time
import os
import numpy as np
from scipy.optimize import curve_fit
import scienceplots

import matplotlib.pyplot as plt
import scienceplots

# plt.style.use(["science", "high-contrast"])
# plt.rcParams["figure.figsize"] = (3.5, 2.6)
# plt.rcParams["font.size"] = "5.4"
# plt.rcParams["axes.labelsize"] = "13"
# plt.rcParams["xtick.labelsize"] = "13"
# plt.rcParams["ytick.labelsize"] = "13"
# plt.rcParams["lines.linewidth"] = "0.9"


from Vacuum_Polarization import Vacuum_Polarization

def linear(x, a):
    return a * (x-1/2)

def parabola(x, a, b,c ):
    return a * x ** 2 + b * x + c

fig1, axRho = plt.subplots()
fig2, axA0 = plt.subplots()
fig3, axmaxA0 = plt.subplots()

pendiente = []
modeArray = []

vp = Vacuum_Polarization(lambdaMin=1., bisectionTol=1e-5, 
        subtractPertModes=True,
        pointsPerMode=8, 
        smoothing=False, 
        hadamard=False,
        )

lambdaValue = vp.lambdaValue

while True:
    try:
        vp.setConfigFromDict(index=0)
    except TypeError as e:
        print(e)
        break
    
    # arr = [10,0,10,0,6,0, 10,0,10]
    # vp.convolveRho(arr)
    # vp.convolveRho(arr)
    p, _ = curve_fit(linear, vp.z, [ float(x) for x in vp.rho])

    vp.rho = vp.rho #- vp.e**2 / np.pi * vp.lambdaValue * (vp.z - 1/2)

    window = vp.nPoints // (vp.maxN + 1)  *2
    kernel = [ np.exp(-x**2/2/window**2) for x in range(-3*window , 3*window+1)]
    # vp.convolveRho(kernel)
    vp.calculateA0Induced()

    # pendiente.append(*p)
    pendiente.append(vp.A0Induced(1))
    
    modeArray.append(vp.maxN)


    axRho.plot(vp.z,vp.rho, '--', label=f"{vp.maxN} modes")
    # axRho.plot(vp.z, linear(vp.z, p) + vp.e**2 / np.pi * vp.A0(vp.z), "--", alpha=0.7)
    axA0.plot(vp.z,vp.A0Induced(vp.z), '--', label=f"{vp.maxN} modes")

    vp.directory = f"Ambjorn N {vp.maxN}"

vp.rho = vp.perturbativeTotalVacuumPolarization(vp.z) - vp.e**2 / np.pi * vp.lambdaValue * (vp.z-1/2)
vp.calculateA0Induced()

axRho.plot(vp.z, vp.rho, label=r"Perturbative $\rho$")
axRho.set_ylabel(r"$\rho$")
axRho.set_xlabel(f"$z$")

axA0.plot(vp.z, vp.A0Induced(vp.z), label=r"Perturbative $A_0^\text{br}$")

axA0.set_ylabel(r"$A_0^\text{br}$")
axA0.set_xlabel(f"$z$")

axmaxA0.plot(modeArray, [vp.A0Induced(1)]*len(modeArray), '--', alpha=0.6,  label=r"Perturbative $A_0^\text{br}$")

# pParabol, _  = curve_fit(parabola, np.array(modeArray), np.array(pendiente))
# axmaxA0.plot(sorted(modeArray), parabola(np.array(sorted(modeArray)), *pParabol), label="curve fit")
axmaxA0.plot(modeArray, pendiente, 'o')
axmaxA0.set_ylabel(r"$A_0^\text{br}(1)$")
axmaxA0.set_xlabel("\# modes")

axA0.legend()
axRho.legend()
axmaxA0.legend()

fig1.suptitle(f"$\lambda = {vp.lambdaValue}, m = {vp.m}$")
fig2.suptitle(f"$\lambda = {vp.lambdaValue}, m = {vp.m}$")
fig3.suptitle(f"$\lambda = {vp.lambdaValue}, m = {vp.m}$")

# fig3.suptitle(f"Comparison to {round(pParabol[0], 5)} x² + {round(pParabol[1], 2)} x + {round(pParabol[2], 2)} ")

date = time.strftime("%Y%m%d-%H%M%S")

if not os.path.exists(date):
    os.makedirs("figures/" + date)

plt.show()
#fig1.savefig(f"figures/{date}/rhoModeComparison.pdf")
#fig2.savefig(f"figures/{date}/A0InducedComparison.pdf")
#fig3.savefig(f"figures/{date}/A0(1) mode comparison.pdf")
