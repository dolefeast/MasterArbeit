import numpy as np
import matplotlib.pyplot as plt


def plot_different_window_filter(total_charge_density, ax_fields, filter_method):
    from math_objects.modify_A0 import modify_A0

    while True:
        plt.draw()
        while True:
            try:
                val = int(input("Input desired window size: "))
                break
            except ValueError:
                print("Input must be an integer")

        first_filter = filter_method(total_charge_density, int(val) // 3 + 1)
        double_filter = filter_method(first_filter, int(val) // 3 + 1)

        double_filter = filter_method(first_filter, int(val) // 3 + 1)
        double_filter = filter_method(first_filter, int(val) // 3 + 1)

        z = np.arange(len(double_filter)) / len(double_filter)
        ax_fields.clear()
        ax_fields.plot(
            z, double_filter, label=f"Moving average filtering. window={val}"
        )

        potential_perturbation = modify_A0(z, double_filter)
        ax_fields.plot(
            z,
            potential_perturbation.sol(z)[0],
            label="Modifying the electric potential",
        )
        plt.pause(0.01)


def plot_each_eigenstate(eigenstate_array, m, lambda_value):
    from math_objects.perturbative_solutions import dirichlet_eigenstate

    fig, ax_array = plt.subplots(len(eigenstate_array))

    for ax, eigenstate in zip(ax_array, eigenstate_array):
        omega_n = eigenstate.p[0]
        ax.plot(eigenstate.x, eigenstate.y[0])
        ax.plot(
            eigenstate.x, dirichlet_eigenstate(eigenstate.x, omega_n, m, lambda_value)
        )
        ax.set_title(f"$\omega_n = {omega_n}$")


def scatter_omegas(eigenvalue_array, ax_omegas, m):
    omega_solution_array = []
    for i, omega_solution in enumerate(eigenvalue_array):
        if omega_solution < 0:
            omega_scatter_color = "r"
            omega_scatter_marker = "x"
        elif omega_solution > 0:
            omega_scatter_color = "g"
            omega_scatter_marker = "o"

        ax_omegas.scatter(
            np.sqrt(omega_solution**2 - m**2) / np.pi + 0.3 * (omega_solution < 0),
            np.abs(omega_solution),
            marker=omega_scatter_marker,
            color=omega_scatter_color,
        )

        omega_solution_array.append(omega_solution)
    return omega_solution_array


def plot_from_0_to_1(array):
    """Given an array [f(x_i)], returns the x such that the plot goes from 0 to 1"""
    x = np.arange(len(array)) / len(array)
    return x, array


def plot_eigenstates(self, modes=(-2, -1, 1, 2), axes=None):
    """

    Parameters:

    """
    N = len(self.eigenstate_array)
    for i, mode_index in enumerate(modes):
        if mode_index < 0:
            shift = -1
        elif mode_index > 0:
            shift = 0
        else:
            print("mode_index=0 is not a valid index")
            continue
        ax[i].plot(
            self.z, self.eigenstate_array[N // 2 + mode_index + shift], alpha=alpha
        )
