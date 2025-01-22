from Vacuum_Polarization import Vacuum_Polarization
import matplotlib.pyplot as plt  

def fixed_point_plot(fun, ax):

    def inner():
        a1 = compute.A0_induced(1)
        fun()
        ax.plot(a1, compute.A0_induced(1), 'o', color=compute.color, label=f"$\lambda={compute.lambda_value}$")

    return inner

compute = Vacuum_Polarization(max_N=1,
        lambda_min=10,
        save=False,
        save_plots=False,
        show_plots=True,
        smoothing=True,
        ambjorn=False,
        m=0,
        relax_parameter=0.6,
        )    

fig, ax = plt.subplots()

compute.core_iteration = fixed_point_plot(compute.core_iteration, ax)

try:
    compute.full_script()
except KeyboardInterrupt:
    pass

ax.set_ylabel("$A_{k+1}$ induced ")
ax.set_xlabel("$A_{k} $ induced ")


handles, labels = plt.gca().get_legend_handles_labels()
by_label = dict(zip(labels, handles))
plt.legend(by_label.values(), by_label.keys())

# plt.show()
