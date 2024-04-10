import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from scripts.plotting import plot_from_0_to_1
from scripts.read_files import read_files

def exp(x, a, b,c):
    return a*np.exp(b*x)+c

max_A0_array = []
 
# Read from the saved solutions
eigenvalue_list,eigenstate_list,  eigenstate_gradient_list, A0_list, rho_list, lambda_list = read_files(m=3.0).values()

max_lambda = max(lambda_list)
print(lambda_list[-2])
fig, (ax_A0, ax_rho) = plt.subplots(2)

for i, (
        A0_induced_file,
        rho_file,
        lambda_value
        ) in enumerate(zip(
                A0_list, 
                rho_list, 
                lambda_list
                )
            ):
    try:
        A0_induced = np.fromfile(
                A0_induced_file, 
                dtype=float,
                sep="\n"
                )
    except Exception:
        continue

    try:
        rho = np.fromfile(
                rho_file, 
                dtype=float, sep="\n"
                )
    except Exception:
        continue

    if i==0:
        A0_min = np.copy(A0_induced)

    alpha = 1- 0.95 + 0.7 * ((lambda_value)/max_lambda)**20
    ax_A0.plot(
            *plot_from_0_to_1(
                A0_induced
                /A0_min
                ), 'b', alpha=alpha), ax_rho.plot(*plot_from_0_to_1(rho),'b',  alpha=alpha)


#ax_A0.plot([], [], 'bo', label="max($A_0(z)$)")

ax_A0.set_title('$A_0^\lambda(z)/A_0^0(z)$')
ax_rho.set_title(r'$\rho(z)$')


ax_A0.legend(loc='best')
plt.show()
