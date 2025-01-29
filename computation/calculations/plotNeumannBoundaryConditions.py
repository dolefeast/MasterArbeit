from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

from Vacuum_Polarization import Vacuum_Polarization

import scienceplots
plt.style.use(["science", "high-contrast"])
omegaArray = np.linspace(-10, 10, 40)

fig, ax = plt.subplots(3, 3, figsize=(10,8))

for m in range(0, 9):
    massIndex = (m//3, m%3)
    for lambdaValue in np.linspace(.1, 5, 5):
        compute = Vacuum_Polarization(m=m, lambdaMin=lambdaValue)

        parametrizedODE = lambda omega: solve_ivp(lambda z, y: compute.KleinGordonEquation(z, y, omega),
                t_span=(0, 1),
                y0=(0, 1), 
                dense_output=True)

        ax[massIndex].plot(omegaArray, [parametrizedODE(omega).sol(1)[0] for omega in omegaArray])
    ax[massIndex].set_title(f"m={compute.m}")
    ax[massIndex].plot(omegaArray, np.zeros_like(omegaArray), "--", alpha=0.5)
    ax[massIndex].plot([-m, m], [0,0], "x", alpha=0.5)


plt.legend(loc="best")
fig.suptitle("$\phi'_\omega(1)$ for different values of m, $\lambda$")
plt.show()
