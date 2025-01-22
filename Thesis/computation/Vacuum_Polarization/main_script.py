from scipy.interpolate import CubicSpline
import numpy as np

def core_iteration(self):
    """
    Calculates and normalizes the eigenstates for the current A0.
    Calculates the vacuum polarization corresponding to this solution family, 
    filters it if necessary, and finally calculates the corresponding induced potential to this A0_induced.
    """
    self.calculate_eigenstates()
    self.normalize_eigenstates()

    rho_new = self.calculate_rho()

    if self.smoothing and not self.ambjorn:
        # rho = filter_rho(rho)

        rho_new += self.e**2/np.pi * self.A0(self.z)
        if self.bcs == "neumann":
            _, rho_derivative_new = self.extend_and_filter(self.z, CubicSpline(self.z, rho_new)(self.z))
        elif self.bcs == "dirichlet":
            _, rho_new = self.extend_and_filter(self.z, rho_new)
        else:
            raise RuntimeError("nothing was done")

    if self.n>70 and not self.n%25 and False:
        print("Averaging..")
        self.rho = (rho_new + self.rho)/2 
    else:
        self.rho = rho_new

    self.A0_induced = self.calculate_A0_induced(self.rho)

def constant_lambda_iterations(self):
    """
    Calculates until convergence or until a certain number of iterations, the backreaction iteration.
    During this calculation it might save different plots for each of the iterations.
    """
    rho_error=1

    if self.save_plots and self.plot_for_each_lambda:
        self.plot_intermediate_steps_setup()

    self.n=0

    self.constant_lambda_A0_list = [ self.A0_history[-1][-1] ]
    self.constant_lambda_A0_induced_list = [ self.A0_induced_history[-1][-1] ]

    while self.n < self.max_n_iterations:

        c = self.calculate_relax_parameter()

        self.A0_induced = self.constant_lambda_A0_induced_list[-1]
        self.A0 = lambda z:  (
                (1-c) * self.constant_lambda_A0_list[-1](z) + c * (- self.lambda_value * (z-1/2) + self.a * self.A0_induced(z)) 
                )

        if isinstance(self.n_iterations, int):
            if self.n >= self.n_iterations: break
        elif self.n_iterations is None:
            if rho_error < self.iterate_tol: 
                break
        else: 
            raise Exception("n_iterations must be either int or None")

        if self.n==0:
            print("inside the while loop")
            print("self.A0(1) = ", self.A0(1))
            print("self.A0_induced(1) = ", self.A0_induced(1))
        # The main calculation. Calculate KG solutions, the associated rho and the associated A0_induced
        self.core_iteration()
        
        # If we want to plot things
        self.plot_intermediate_steps() 

        if self.n_iterations is None:
            try:
                # Tries to calculate the error if max_prev_rho is defined
                rho_error = abs(max(self.rho) - max_prev_rho) / abs(max_prev_rho)
                max_prev_rho = max(self.rho)
            except NameError:
                # If it wasn't, define it and calculate again
                max_prev_rho = max(self.rho)

        self.A0 = CubicSpline(self.z, self.A0(self.z))

        self.constant_lambda_A0_list.append(self.A0)
        self.constant_lambda_A0_induced_list.append(self.A0_induced)
        self.n += 1
    else:
        if self.save_plots:
            self.fig4.savefig(
                    "figures/"+self.directory+"/"+self.bcs+"/intermediate_steps_lambda_value_"+self.float_to_str(self.lambda_value)+".png"
                    )
        raise RecursionError("max_n_iterations reached")

    if self.save_plots: 
        self.fig4.savefig(
                "figures/"+self.directory+"/"+self.bcs+"/intermediate_steps_lambda_value_"+self.float_to_str(self.lambda_value)+".png")

def walkback(self):
    """
    If a failure in the calculation was found, "walk back" i.e. return to the last lambda_value that lead to converged solutions, and reduce 1. the lambda_step 2. the relax parameter in the update law
    """

    self.lambda_value -= self.lambda_step 
    self.lambda_step *= self.walkback_factor
    self.relax_parameter *= self.walkback_factor
    self.lambda_value += self.lambda_step

def full_script(self):

    while self.lambda_value < self.lambda_max:
        if self.lambda_step < self.lambda_step_min:
            break

        print(15*"#")
        print('λ =', self.lambda_value)
        self.color = next(self.color_cycle) # To distinguish between different lambda values

        try:
            self.constant_lambda_iterations() 

            # If no exception was raised, a self-consistent solution was found and it converged
            if self.plot:
                self.eigenvalues_array.append(self.eigenvalues)
                self.lambda_array.append(self.lambda_value)
                self.plot_rho_A0_induced(self.rho, self.A0_induced(self.z), '', color=self.color, alpha=1)

            print(f'Converged in n={self.n} iterations. A0_induced(1) =', self.A0_induced(1))
            if self.save:
                self.save_solutions()

        except RecursionError as ex:
            break
        except RuntimeError as exception:
            self.walkback()

            print(exception) # To know if it didn't converge or if no solution was found
            print(f'Error found at n={self.n}, A0_induced(1) =', self.A0_induced(1))

            continue
        except Exception as exception:
            break
        except KeyboardInterrupt as exception:
            break

        self.A0_history.append(self.constant_lambda_A0_list)
        self.A0_induced_history.append(self.constant_lambda_A0_induced_list)


        self.lambda_value += self.lambda_step

    
    if self.save_plots:
        self.save_all_plots()
    if self.show_plots:
        self.ax_eigenvalues.plot(self.lambda_array, self.eigenvalues_array, 'b')
        self.plt.show()
