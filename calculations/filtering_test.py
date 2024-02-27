import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scripts.filtering import moving_average
from scripts.plotting import plot_from_0_to_1

unfiltered_charge_density = np.genfromtxt('./saved_solutions/lambda_1_mass_1.csv', delimiter=",")

x, y = unfiltered_charge_density

# To fix the fact that the information of where the original boundaries are, create x with len(x) = len(y), but ranging from the corresponding places
# example: y has  len=500.
# if the periodic y is done extending it by half an oscillation in each direction, it now has len = 500 + 250 + 250=1000.
# original x ranged from 0 to 1.
# now it has to range from -0.5 to 1.5.
# only keep values of y for which x is in betweeen 0 and 1. 

y_list = list(y)
len_y_list = len(y_list)
periodic_charge_density = y_list[3*len_y_list//4:-1] + y_list + y_list[:len_y_list//4]

convolution = moving_average(periodic_charge_density, 450)
n_points = len(convolution)

convolution = convolution[n_points//5:4*n_points//5]
x, y = plot_from_0_to_1(convolution)
convolution = sp.interpolate.UnivariateSpline(x, y)
# y = convolution(x)
zeros = convolution.roots()

idx = np.where(
            np.logical_and(
                x > min(zeros),
                x < max(zeros)
            )
            )
x = x[idx]
y = y[idx]
# plt.plot(*plot_from_0_to_1(periodic_charge_density))
# plt.plot(*plot_from_0_to_1(convolution))
plt.plot(x, y)
plt.plot(zeros, 3*[0], 'o')


plt.show()
