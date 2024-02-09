import matplotlib.pyplot as plt
from math_objects.savitzky_golay import savitzky_golay
import numpy as np
import scipy as sp

n_points = 500
x = np.linspace(0.01, 0.99, n_points)
y = x*(1-x)/np.tan(np.pi*x)

N = 1000
series = 0
for n in range(-N, N):
    series = series + np.sin(np.pi*abs(n)*x)*np.cos(np.pi*x*n)
    # series = series + savitzky_golay(np.sin(np.pi*n*x)*np.cos(np.pi*x*n), 40, 1)
    cut_off_series = series*x*(1-x)
    # plt.plot(x, cut_off_series)

plt.plot(x, y/2, label="analytical solution") # a factor of four is picked up in the series
plt.plot(x, cut_off_series, label="Cut off series (unfiltered)", alpha=0.7, linewidth=0.7)

filtered = savitzky_golay(cut_off_series, int(n_points/1000 * 90 )+ 1, 3)
plt.plot(x, filtered, label="Cut off series (filtered)")

plt.legend(loc="best")
plt.show()
