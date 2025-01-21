import matplotlib.pyplot as plt

def plotting_setup(self):
    self.fig1 = self.plt.figure(r"Vacuum polarization")
    self.fig2 = self.plt.figure(r"A_0 induced")
    self.fig3 = self.plt.figure("Mode energy evolution")
    self.fig4 = self.plt.figure("Values at intermediate steps")

    self.ax_rho = self.fig1.subplots(1)
    self.ax_A0_induced = self.fig2.subplots(1)
    self.ax_eigenvalues = self.fig3.subplots(1)
    self.ax_relax, self.ax_intermediate_omegas = self.fig4.subplots(2)

    self.eigenvalues_array = []
    self.lambda_array = []

def plot_rho_A0_induced(self, rho, A0_induced, fmt, color, alpha,):
    self.ax_rho.plot(self.z, rho, fmt, color=self.color, alpha=alpha)
    self.ax_A0_induced.plot(self.z, A0_induced, fmt, color=self.color, alpha=alpha)

def plot_intermediate_steps_setup(self):
    self.fig4.suptitle(r"$\lambda$={}".format(round(self.lambda_value, 4)))

    self.ax_relax.cla()
    self.ax_intermediate_omegas.cla()

    self.ax_intermediate_omegas.set_xlabel(r"n")
    self.ax_relax.set_ylabel(r"$A_0(1) +$ {}".format(round(self.lambda_value/2, 5)))
    self.ax_intermediate_omegas.set_ylabel(r"$\omega$")

def plot_intermediate_steps(self):
    if self.plot:
        self.plot_rho_A0_induced(self.rho, self.A0_induced(self.z), '-', color=self.color, alpha=0.2)

    if self.plot_for_each_lambda and self.save_plots:
        self.ax_relax.plot(self.n, self.A0(1) + round(self.lambda_value/2, 2), 'x', color=self.color)
        self.ax_intermediate_omegas.plot(self.n, self.eigenvalues[self.max_N], 'x', color=self.color)


def save_all_plots(self):
    self.ax_eigenvalues.plot(self.lambda_array, self.eigenvalues_array, 'b')
    self.fig1.savefig("figures/"+self.directory+"/"+self.bcs+"/vacuum_polarization.png")
    self.fig2.savefig("figures/"+self.directory+"/"+self.bcs+"/A0_induced.png")
    self.fig3.savefig("figures/"+self.directory+"/"+self.bcs+"/eigenvalues.png")
