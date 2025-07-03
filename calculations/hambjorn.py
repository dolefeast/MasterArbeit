import matplotlib.pyplot as plt
import scienceplots
import numpy as np

from Vacuum_Polarization import Vacuum_Polarization

plt.style.use(["science", "high-contrast"])
plt.rcParams["figure.figsize"] = (3.9, 2.9)
plt.rcParams["font.size"] = "5.4"
plt.rcParams["axes.labelsize"] = "13"
plt.rcParams["xtick.labelsize"] = "13"
plt.rcParams["ytick.labelsize"] = "13"
plt.rcParams["lines.linewidth"] = "0.9"

vp = Vacuum_Polarization(lambdaMin = 1)

fig, ax = plt.subplots()

plt.plot(vp.z, vp.perturbativeTotalVacuumPolarization(vp.z) - vp.e**2 / np.pi * vp.lambdaValue *(vp.z-1/2), label="With parallel transport")
plt.plot(vp.z, vp.perturbativeTotalVacuumPolarization(vp.z), label="Without parallel transport")
plt.plot(vp.z, vp.perturbativeModeRho(1,vp.z), label=r"$\frac{1}{2}(\rho_1(z) + \rho_{-1}(z))$")


ax.legend()
ax.set_ylabel(r"$\rho$")
ax.set_xlabel(r"$z$")

fig.tight_layout(rect=[0, 0.03, 1, 0.95])
fig.suptitle(f"m = {vp.m}, $\lambda$ = {vp.lambdaValue}")

 
fig.savefig("hambjorn.pdf")