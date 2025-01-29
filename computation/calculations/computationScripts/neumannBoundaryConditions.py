from Vacuum_Polarization import Vacuum_Polarization
import matplotlib.pyplot as plt
import numpy as np

compute = Vacuum_Polarization(max_N=40,
        lambda_min = 1,
        smoothing=False,
        bcs="neumann"
        )

compute.n = 0

# compute.core_iteration()
# np.savetxt( "neumann_vacuum_pol.txt", compute.rho)
# plt.plot(compute.rho)
# plt.show()
# exit()

rho = []
with open("neumann_vacuum_pol.txt") as rho_text:
    for i in rho_text:
        rho.append(float(i))

rho = np.array(rho)
rho_old = np.copy(rho)

rho = compute.filter_rho(rho)

plt.plot(rho_old, label=r"$\rho$ old")
plt.plot(rho, label=r"$\rho$ filtered")

plt.legend(loc="best")
plt.show()
