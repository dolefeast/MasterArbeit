import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from pathlib import Path
from scripts.plotting import plot_from_0_to_1

things = ['A0_field', 'rho']
A0_list, rho_list = [
        list(Path(f'saved_solutions/dirichlet/{thing}').glob('*lambda_*_mass_4_0*'))
        for thing in things
        ]

A0_list = sorted(A0_list)[::-1]
rho_list = sorted(rho_list)[::-1]


fig, (ax_A0, ax_rho) = plt.subplots(2)
for i, (rho_file, A0_perturbation_file) in enumerate(zip(A0_list, rho_list)):
    try:
        A0_perturbation = np.fromfile(A0_perturbation_file, dtype=float, sep='\n')
    except Exception:
        continue

    try: 
        rho = np.fromfile(rho_file, dtype=float, sep='\n')
    except Exception:
        continue

    alpha = 0.9 - 0.8 * ((i+1)/len(A0_list))**1
    ax_A0.plot(*plot_from_0_to_1(A0_perturbation), 'b', alpha=alpha)
    ax_rho.plot(*plot_from_0_to_1(rho),'b',  alpha=alpha)


ax_rho.set_title('$A_0(z)$')
ax_A0.set_title(r'$\rho(z)$')
plt.show()
