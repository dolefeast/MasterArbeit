import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import regex as re

from pathlib import Path
from scripts.plotting import plot_from_0_to_1

from scripts.float_to_str import str_to_float

things = ['A0_field', 'rho']
A0_list, rho_list = [
        list(Path(f'saved_solutions/dirichlet/{thing}').glob('*lambda_*_mass_3_0*'))
        for thing in things
        ]

lambda_re = re.compile("lambda_\d+_\d+")
float_re = re.compile("\d+_\d+")
lambda_list = [str_to_float(float_re.findall(str(filename))[0]) for filename in A0_list]
#mass_list = [float_re.findall(str(filename))[0] for filename in A0_list]

max_lambda = max(lambda_list)
fig, (ax_A0, ax_rho) = plt.subplots(2)
for i, (rho_file, A0_perturbation_file, lambda_value) in enumerate(zip(A0_list, rho_list, lambda_list)):
    try:
        A0_perturbation = np.fromfile(A0_perturbation_file, dtype=float, sep='\n')
    except Exception:
        continue

    try: 
        rho = np.fromfile(rho_file, dtype=float, sep='\n')
    except Exception:
        continue

    alpha = 1- 0.9 + 0.8 * ((lambda_value)/max_lambda)**20
    ax_A0.plot(*plot_from_0_to_1(A0_perturbation/lambda_value), 'b', alpha=alpha)
    ax_rho.plot(*plot_from_0_to_1(rho),'b',  alpha=alpha)


ax_rho.set_title('$A_0(z)$')
ax_A0.set_title(r'$\rho(z)/\lambda$')
plt.show()
