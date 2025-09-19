import matplotlib.pyplot as plt
import numpy as np

from Vacuum_Polarization import Vacuum_Polarization

import matplotlib_format

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

vp = Vacuum_Polarization(hadamard=True)

vp.setConfigFromDict()

ax1.plot(vp.z, vp.rho, label=f"{vp.maxN} modes") # [vp.nPoints:2*vp.nPoints + 1])
ax2.plot(vp.z, vp.A0Induced(vp.z), label=f"{vp.maxN} modes") # [vp.nPoints:2*vp.nPoints + 1])

plt.show()