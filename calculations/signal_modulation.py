import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from scripts.filtering import moving_average
from scripts.bao_filtering import bao_filtering

from scripts.plotting import plot_from_0_to_1





def main():
    unfiltered_charge_density = np.genfromtxt('./saved_solutions/lambda_1_mass_1.csv', delimiter=",")

    x, y = unfiltered_charge_density 

    plt.plot(x, y, alpha= 0.4)
    x,y = bao_filtering(x, y, size=0.06)
    plt.plot(x, y)

    # plt.plot(x, no_wiggle_y(x))

    plt.show()

if __name__ == "__main__":
    main()
