import matplotlib.pyplot as plt

from snippets.compare_renormalization_fields_induced import main
from plot.format_ax_rho_A0_induced import format_ax_rho_A0_induced

# fig, (ax_rho, ax_A0_induced) = plt.subplots(2)

# format_ax_rho_A0_induced(ax_rho, ax_A0_induced)

main(
        m=0,
        a = 1,
        # ax_rho=ax_rho,
        # ax_A0_induced=ax_A0_induced,
        lambda_value = 5,
        directory="a_evolution",
        max_alpha = 1
        )

# plt.show()
