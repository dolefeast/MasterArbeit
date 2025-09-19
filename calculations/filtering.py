import matplotlib.pyplot as plt
import numpy as np

from Vacuum_Polarization import Vacuum_Polarization
import matplotlib_format

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

vp = Vacuum_Polarization(hadamard=True)
while True:
    # Reading data
    try:
        vp.setConfigFromDict()

        vp.hadamard = False
        vp.rho = vp.rho - vp.lambdaValue*(vp.z-1/2)/np.pi

        # vp.rho = np.concatenate((vp.rho, vp.rho, vp.rho ))

        # plt.plot(vp.rho/10, label="before filtering") #[vp.nPoints:2*vp.nPoints + 1])

        # Two oscillations
        window = vp.nPoints // (vp.maxN + 1) *2

        kernel = [ np.exp(-x**2*4/window**2)/7 for x in range(-3*window , 3*window+1)]

        vp.convolveRho()

        ax1.plot(vp.z, vp.rho, label=f"{vp.maxN} modes") # [vp.nPoints:2*vp.nPoints + 1])
        ax2.plot(vp.z, vp.rho, label=f"{vp.maxN} modes") # [vp.nPoints:2*vp.nPoints + 1])
        ax2.plot(vp.z[:len(kernel)], kernel, label="Convolution kernel") # [vp.nPoints:2*vp.nPoints + 1])
    except TypeError:
        break


    # plt.plot(vp.perturbativeTotalVacuumPolarization(vp.z) - vp.lambdaValue * (vp.z- 1/2)/np.pi, label="perturbative")

ax2.set_xlabel(r'$z$')
ax2.set_xlim((0,7*window/vp.nPoints))

ax1.plot(vp.z, vp.perturbativeTotalVacuumPolarization(vp.z) - vp.lambdaValue/np.pi * (vp.z - 1/2), "--", label=r"Perturbative $\rho$")

ax1.set_xlabel(r'$z$')
ax1.set_ylabel(r'$\rho$')

ax2.legend()
ax1.legend()


fig1.tight_layout()
fig2.tight_layout()

fig1.savefig("figures/unfilteredRho.pdf")
fig2.savefig("figures/unfilteredRho+kernel.pdf")


