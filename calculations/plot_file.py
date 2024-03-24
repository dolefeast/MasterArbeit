import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from pathlib import Path
from scripts.plotting import plot_from_0_to_1

things = ['A0_field', 'charge_density']
A0_list, charge_density_list = [
        list(Path(f'saved_solutions/dirichlet/{thing}').glob('*lambda_*_mass_1_12*'))
        for thing in things
        ]

fig, (ax_A0, ax_charge_density) = plt.subplots(2)
for i, (charge_density_file, A0_perturbation_file) in enumerate(zip(A0_list, charge_density_list)):
    try:
        A0_perturbation = np.fromfile(A0_perturbation_file, dtype=float, sep='\n')
    except Exception:
        continue

    try: 
        charge_density = np.fromfile(charge_density_file, dtype=float, sep='\n')
    except Exception:
        continue

    alpha = 0.2 + 0.8 * ((i+1)/len(A0_list))**4
    ax_A0.plot(*plot_from_0_to_1(A0_perturbation), alpha=alpha)
    ax_charge_density.plot(*plot_from_0_to_1(charge_density), alpha=alpha)


ax_charge_density.set_title('$A_0(z)$')
ax_A0.set_title(r'$\rho(z)$')
plt.show()
