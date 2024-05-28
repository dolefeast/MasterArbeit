import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots()

data = np.genfromtxt('lambda_c_convergence2.txt', delimiter=',')


ax.plot(*data,
        'o',
        )

plt.show()
