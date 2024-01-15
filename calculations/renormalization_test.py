import matplotlib.pyplot as plt
from math_objects.savitzky_golay import savitzky_golay
import numpy as np
import scipy as sp

n_points = 500
x = np.linspace(0.01, 0.99, n_points)
y = x*(1-x)/np.tan(np.pi*x)

N = 1000
series = 0
for n in range(N):
    series = series + savitzky_golay(np.sin(np.pi*n*x)*np.cos(np.pi*x*n), 40, 1)

plt.plot(x, y/4, label="analytical solution") # a factor of for is picked up in the series
cut_off_series = series*x*(1-x)
plt.plot(x, cut_off_series, label="Cut off series (unfiltered)")

filtered = savitzky_golay(cut_off_series, 40, 1)
plt.plot(x, filtered, label="Cut off series (filtered)")

plt.legend(loc="best")
plt.show()
