import scipy as sp
import matplotlib.pyplot as plt
import numpy as np


# y'' + y = 0
# y(0) = y(1) = 0

n_points = 50


def ODE(x, y):
    ode = np.vstack((y[1], -100 * y[0] + 10 * y[1] + 5))
    print(ode.shape)
    return ode


def boundary_conditions(ya, yb):
    return np.array((ya[1] - ya[0], yb[1]))


x = np.linspace(0, 1, n_points)

y_a = np.ones((2, x.size))
print(np.shape(y_a))
# y_b = np.zeros((2, x.size))

res = sp.integrate.solve_bvp(ODE, boundary_conditions, x, y_a)

plt.plot(x, res.sol(x)[0], label="solution")
plt.plot(x, res.sol(x)[1], label="derivative")

plt.legend(loc="best")
plt.show()
