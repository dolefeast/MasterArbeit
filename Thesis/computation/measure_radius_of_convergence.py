from Vacuum_Polarization import Vacuum_Polarization
from scipy.interpolate import CubicSpline
import numpy as np
import matplotlib.pyplot as plt  

import os
import sys

def get_ax_lims(ax):
    return (ax.get_xlim(), ax.get_ylim())

def fixed_point_plot(fun):

    def inner():

        A0_induced_value = compute.A0_induced(1)
        A0_full_value = compute.lambda_value/2 + compute.A0(1)
        fun()
        fixed_point_line_induced.append((A0_induced_value, compute.A0_induced(1)))
        fixed_point_line_full.append((A0_full_value, compute.lambda_value/2 + compute.A0(1)))
        # ax.plot(a1, compute.A0_induced(1), 'o', color=compute.color, label=f"$A0_induced(1)={compute.A0_induced(1)}$")

    return inner

compute = Vacuum_Polarization(
        max_N=int(sys.argv[1]),
        iterate_tol = 0.0005,
        lambda_min=10,
        lambda_max=17,
        lambda_step=6,
        save=False,
        save_plots=False,
        show_plots=False,
        smoothing=True,
        ambjorn=False,
        m=0,
        relax_parameter=0.6,
        max_n_iterations=500,
        )    

# I know lambda_value = 16 is a stable solution
# want to measure the radius of convergence
# So: Get to lambda_value = 16, store the configuration
# of the system in that case, and then study what happens
# when using different starting A0_induced

try:
    compute.full_script()
except KeyboardInterrupt:
    pass

# It should have stopped at a configuration that yields a bifurcation
# because it breaks in case of RecursionError

bifurcation_config = dict()
bifurcation_config["eigenvalues"] = compute.eigenvalues
bifurcation_config["eigenstates"] = compute.eigenstates
bifurcation_config["A0_induced"] = compute.A0_induced(compute.z)

bifurcation_config["lambda_value"] = 17

# Looking for the radius of convergence of a certain configuration
# This means: What's the furthest away I can start from the converged
# A_0 so that my problem converges?

relax_parameter_array = [0.001, 0.01, 0.05, 0.1][::-1]
A0_factor_array = np.arange(-15, 15, 2)

for c in relax_parameter_array:
    print(f"c = {c}")

    success_array = []
    fig1, (ax_A0_induced_fixed_point, ax_A0_induced_iteration) = plt.subplots(nrows=1, ncols=2, figsize=(16,9))
    fig2, (ax_A0_full_fixed_point, ax_A0_full_iteration) = plt.subplots(nrows=1, ncols=2, figsize=(16,9))
    fig1.suptitle(f"c = {c}")
    fig2.suptitle(f"c = {c}")

    compute.relax_parameter = c
    compute.iterate_tol = c/10
    compute.core_iteration = fixed_point_plot(compute.core_iteration)

    for A0_factor in A0_factor_array:
        directory = f"convergence_radius_max_N_{compute.max_N}"

        fixed_point_line_induced = []
        fixed_point_line_full = []

        compute.set_config_from_dict(bifurcation_config)
        compute.A0_induced_history = [
                [
                CubicSpline(
                compute.z, 
                A0_factor * bifurcation_config["A0_induced"]
                )
                    ]
                ] # Need to overwrite this

        try:
            compute.color = next(compute.color_cycle)
            compute.constant_lambda_iterations()
            success = compute.A0_induced(1)
        except KeyboardInterrupt:
            break
        except RuntimeError as e: 

            print(e)
            success = False
        finally:
            if len(fixed_point_line_induced) != 0:
                ax_A0_induced_iteration.plot([x for x, y in fixed_point_line_induced], 
                        'x-',
                        label=f"A0 induced(1) = {round(A0_factor * bifurcation_config['A0_induced'][-1], 3)}",
                        color=compute.color)

                ax_A0_induced_fixed_point.plot(
                        [x  for x, y in fixed_point_line_induced],
                        [y  for x, y in fixed_point_line_induced], 
                        'x-',
                        label=f"A0 induced(1) = {round(A0_factor * bifurcation_config['A0_induced'][-1], 3)}",
                        color=compute.color)

                ax_A0_full_iteration.plot([x for x, y in fixed_point_line_full], 
                        'x-',
                        label=f"A0 full(1) = {round(fixed_point_line_full[0][0], 3)}",
                        color=compute.color)

                ax_A0_full_fixed_point.plot(
                        [x  for x, y in fixed_point_line_full],
                        [y  for x, y in fixed_point_line_full], 
                        'x-',
                        label=f"A0 full(1) = {round(fixed_point_line_full[0][0], 3)}",
                        color=compute.color)

        
        success_array.append(success)
    ax_A0_induced_iteration.set_xlabel("k")
    ax_A0_induced_iteration.set_ylabel("$A_k(1)$ induced")

    ax_A0_full_iteration.set_xlabel("k")
    ax_A0_full_iteration.set_ylabel(f"full $A_k(1)$ + {compute.lambda_value/2}")

    ax_A0_induced_fixed_point.set_xlabel("$A_k(1)$ induced")
    ax_A0_induced_fixed_point.set_ylabel("$A_{k+1}(1)$ induced")

    ax_A0_full_fixed_point.set_xlabel("full $A_{{k}}(1)$ + {}".format(compute.lambda_value/2))
    ax_A0_full_fixed_point.set_ylabel("full $A_{{k+1}}(1)$ + {}".format(compute.lambda_value/2))

    ax_A0_induced_iteration.legend(loc="best")
    ax_A0_full_iteration.legend(loc="best")
    ax_A0_induced_fixed_point.legend(loc="best")
    ax_A0_full_fixed_point.legend(loc="best")

    ax_A0_induced_iteration.legend(loc="best")
    ax_A0_full_iteration.legend(loc="best")
    ax_A0_induced_fixed_point.legend(loc="best")
    ax_A0_full_fixed_point.legend(loc="best")

    ax_fixed_point_list = [ax_A0_induced_fixed_point, ax_A0_full_fixed_point]

    for ax in ax_fixed_point_list:
        ax_lims = get_ax_lims(ax) 
        print(ax_lims)
        ax.plot([-999, -999], [999, 999], 'g--', alpha=0.2)
        ax.set_xlim(ax_lims[0])
        ax.set_ylim(ax_lims[1])

    try:
        os.mkdir(f"figures/{directory}")
    except FileExistsError:
        pass

    fig1.savefig(f"figures/{directory}/fixed_point_A0_induced_c_{compute.float_to_str(c)}.pdf")
    fig2.savefig(f"figures/{directory}/fixed_point_A0_full_c_{compute.float_to_str(c)}.pdf")
