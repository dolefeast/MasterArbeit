import numpy as np
from tests.increase_lambda import main
from utils.float_to_str import float_to_str

m = 0
lambda_min = 14.5
# lambda_step = 0.01
lambda_max = 30
n_iterations = None
verbose = 0
plot = False
save_solutions = False
read_solutions = True
tol=5e-3
a=1

lambda_step_inverse_array = np.linspace(50, 100, 15)
lambda_critical_array = []

with open(f'saved_solutions/dirichlet/lambda_convergence_a_{float_to_str(a, 3)}.txt', "a") as results_file:
    for lambda_step_inverse in lambda_step_inverse_array:
        print(f"lambda_step = ({lambda_step_inverse})**-1 = {lambda_step_inverse**-1}") 
        lambda_critical = main(
            m=m,
            a=a,
            lambda_min=lambda_min,
            lambda_step=lambda_step_inverse**-1,
            lambda_max=lambda_max,
            n_iterations=n_iterations,
            verbose=verbose,
            plot=plot,
            save_solutions=save_solutions,
            read_solutions=read_solutions,
            tol=tol,
            directory="a_evolution",
                )
        lambda_critical_array.append(lambda_critical)
        results_file.write(f"{lambda_step_inverse}, {lambda_critical}\n")
