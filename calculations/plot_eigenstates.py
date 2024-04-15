import matplotlib.pyplot as plt

from scripts.read_files import read_files
from scripts.plotting import plot_from_0_to_1

m=2.
data = read_files(m=m)

eigenstate_array_array = data["eigenstate_array"]
eigenvalue_array_array = data["eigenvalue_array"]
A0_induced_array = data["A0_induced"]
lambda_value_array = data["lambda_value"]

max_lambda = max(lambda_value_array)

fig, (ax, ax_anti, ax_total) = plt.subplots(3, figsize=(16, 9))

mode = 0
for eigenstate_array, eigenvalue_array, A0_induced, lambda_value in zip(
        eigenstate_array_array,
        eigenvalue_array_array,
        A0_induced_array,
        lambda_value_array,
        ):
    z, eigenstate = plot_from_0_to_1(eigenstate_array[50 + mode])
    eigenvalue = eigenvalue_array[50 + mode]

    anti_eigenstate = eigenstate_array[50 - mode - 1]
    anti_eigenvalue = eigenvalue_array[50 - mode - 1]

    A0_field = - lambda_value * (z - 1/2) + A0_induced

    alpha = 0.2 + 0.7 * (lambda_value / max_lambda)**4

    rho_1 = (eigenvalue - A0_field) * eigenstate ** 2
    anti_rho_1 = (anti_eigenvalue - A0_field) * anti_eigenstate ** 2
    ax.plot(z,rho_1, 'g', alpha=alpha)
    ax_anti.plot(z, anti_rho_1, 'b', alpha=alpha)
    ax_total.plot(z, rho_1 + anti_rho_1, 'k', alpha=alpha)

fig.suptitle(f"$max\, \lambda = {max_lambda}, m={m}$")  
ax.set_title(r'$(\omega_1-eA_0^{(\lambda)} (z))| \phi_1 (z) |^2$')
ax_anti.set_title(r'$(\omega_{-1}-eA_0^{(\lambda)} (z))| \phi_{-1} (z) |^2$')
ax_total.set_title(r'$\rho_1(z) + \rho_{-1} (z)$')

ax.set_xlabel(r'z')
ax_anti.set_xlabel(r'z')
ax_total.set_xlabel(r'z')
fig.tight_layout()
plt.show()

