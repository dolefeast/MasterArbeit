import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


from scripts.calculate_charge import calculate_charge
from scripts.read_files import read_files

def exp(x, a, b,c):
    return a*np.exp(b*x)+c

m = 3

max_A0_array = []
induced_charge_list = []
 
# Read from the saved solutions
(
        eigenvalue_list,
        a,
        b,
        A0_list,
        rho_list,
        lambda_list
) = read_files(m=m).values()

max_lambda = max(lambda_list)
fig, (ax_A0, ax_q) = plt.subplots(2)

for i, (
        A0_perturbation,
        rho,
        lambda_value
        ) in enumerate(
                zip(
                A0_list, 
                rho_list, 
                lambda_list
                )
            ):

    max_A0_array.append(max(A0_perturbation))

    q = calculate_charge(rho)
    induced_charge_list.append(q)


lambda_list = np.array(lambda_list)

# Since it appears to have exponential dependence for λ>18
idx = np.arange(len(lambda_list)) 
idx = np.where(lambda_list>=18)

# Fitting the data to an exponential
popt, pcov = curve_fit(exp, lambda_list[idx], np.array(max_A0_array)[idx], maxfev=5000)
print(f'Optimal parameters: {list(popt)} for m={m}')

ax_A0.plot(
        lambda_list, 
        max_A0_array, 
        'bo', 
        label="max($A_0(z)$)",
        ) # The max(A0) data
ax_A0.plot(
        lambda_list, 
        exp(lambda_list, *popt),
        label='fit to curve',
        ) # The predicted behaviour

ax_q.plot(
        lambda_list, 
        np.array(lambda_list) - np.array(induced_charge_list),
        ) # The predicted behaviour

ax_A0.set_ylabel('max($A_0(z)$)')
ax_A0.set_xlabel('$\lambda$')

ax_q.set_ylabel(r'$\lambda/e- \int_0^{1/2} \rho dz$')
ax_q.set_xlabel('$\lambda$')
ax_A0.legend(loc='best')
fig.suptitle(r'$m={}$'.format(m))
plt.show()
