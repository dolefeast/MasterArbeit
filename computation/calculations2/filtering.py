import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import numpy as np
from scipy.signal import butter, filtfilt

from numpy.fft import rfft

from Vacuum_Polarization import Vacuum_Polarization

def butter_lowpass_filter(data, cutoff, fs, order):
    normal_cutoff = cutoff
    
    b, a = butter(order, normal_cutoff, btype='low', analog=True)
    y = filtfilt(b, a, data)
    return y

# Reading data
vp = Vacuum_Polarization(hadamard=True)
vp.setConfigFromDict()

vp.hadamard = False

vp.rho = vp.rho #- vp.lambdaValue*(vp.z-1/2)/np.pi

# vp.rho = np.concatenate((vp.rho, vp.rho, vp.rho ))

# plt.plot(vp.rho/10, label="before filtering") #[vp.nPoints:2*vp.nPoints + 1])

# Two oscillations
window = vp.nPoints // (vp.maxN + 1) *2

# kernel = [ np.exp(-x**2*4/window**2) for x in range(-3*window , 3*window+1)]
vp.convolveRho()

plt.plot(vp.rho, label="filtered") # [vp.nPoints:2*vp.nPoints + 1])
plt.plot(vp.perturbativeTotalVacuumPolarization(vp.z) - vp.lambdaValue * (vp.z- 1/2)/np.pi, label="perturbative")

plt.legend()


plt.show()
