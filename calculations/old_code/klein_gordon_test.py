eigenstates = []
omega_solution_array = []

eigenstate_array = system.calculate_N_eigenstates(
    -omega_boundary,
    omega_boundary,
    n_of_solutions,
)


total_charge_density_callable = system.calculate_total_charge_density(eigenstate_array)

for i, eigenstate in enumerate(eigenstate_array):
    if not eigenstate.success:
        fmt = "r"
    else:
        fmt = "g"

    omega_solution = eigenstate.p[0]
    if float_in_array(omega_solution, omega_solution_array):
        continue

    # ax_fields.plot(system.z,  system.calculate_charge_density(eigenstate)(system.z), fmt, alpha=0.5, linewidth=0.6)
    ax_omegas.scatter(i, omega_solution, color=fmt)
    ax_omegas.set_title(f"$\lambda={lambda_value}, m={scalar_mass}$")
    omega_solution_array.append(omega_solution)

total_charge_density_array = total_charge_density_callable(system.z)
total_charge_density_array[0] = 0

ax_fields.plot(
    system.z, total_charge_density_array, label="Rough charge density", alpha=0.6
)

# SMOOTHING
smoothed_charge_density = savitzky_golay(
    total_charge_density_array, 10 * system.n_points // (n_of_solutions + 1), 2
)

ax_fields.plot(
    system.z, smoothed_charge_density, label="Savitzky golay smoothing 1st time"
)

smoothed_charge_density = savitzky_golay(
    smoothed_charge_density, 10 * system.n_points // (n_of_solutions + 1), 4
)
ax_fields.plot(
    system.z, smoothed_charge_density, label="Savitzky golay smoothing 2nd time"
)

# electric_field = system.new_electric_field(eigenstate_array)
vector_field = system.new_vector_field(smoothed_charge_density)
# ax_fields.plot(system.z, vector_field.y[0], label='$A_0$ field')
ax_fields.plot(system.z, vector_field.y[0], label="$A_0$ field")

ax_fields.legend(loc="best")

plt.tight_layout()
plt.show()
