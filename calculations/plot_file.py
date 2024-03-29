import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from scripts.plotting import plot_from_0_to_1
import regex as re

things = ["A0_field", "rho"]
A0_list, rho_list = [
    list(Path(f"saved_solutions/dirichlet/{thing}").glob("*lambda_*_mass_3_0*"))
    for thing in things
]

p = re.compile("lambda_\d+_\d+")
lambda_values = [re.findall(p, str(filename)) for filename in A0_list]
value = re.compile("\d+_\d+")

for i, filename in enumerate(lambda_values):
    lambda_value = re.findall(value, str(filename))[0] 
    lambda_value = float(lambda_value.replace("_", "."))
    lambda_values[i] = lambda_value

max_lambda = max(lambda_values)
fig, (ax_A0, ax_rho) = plt.subplots(2)
for i, (rho_file, A0_perturbation_file, lambda_value) in enumerate(zip(A0_list, rho_list, lambda_values)):
    try:
        A0_perturbation = np.fromfile(A0_perturbation_file, dtype=float, sep="\n")
    except Exception:
        continue

    try:
        rho = np.fromfile(rho_file, dtype=float, sep="\n")
    except Exception:
        continue

    alpha = 1- 0.9 + 0.8 * ((lambda_value) / max_lambda) ** 150
    ax_A0.plot(*plot_from_0_to_1(A0_perturbation), "b", alpha=alpha)
    ax_rho.plot(*plot_from_0_to_1(rho), "b", alpha=alpha)


ax_rho.set_title("$A_0(z)$")
ax_A0.set_title(r"$\rho(z)$")
plt.show()
