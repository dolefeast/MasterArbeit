import matplotlib.pyplot as plt
from scripts.read_files import read_files


masses = [1, 2, 3, 3.5, 4, 5]


for m in masses[1:]:
    data = read_files(m=m)

    lambda_value = data["lambda_value"][-1]
    omega_1 = data["eigenvalue_array"][-1][50]
    A0_induced = data["A0_induced"][-1]

    A0_total = - lambda_value/2 + A0_induced

    plt.plot(
            m, 
            - (omega_1 - A0_total[-1])**2 + m**2,
            'bo'
            )
plt.show()

