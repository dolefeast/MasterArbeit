from tests.lambda_c_convergence import main

m = 0
lambda_min = 16

inverse_lambda_step_min = 80
inverse_lambda_step_max = 120
inverse_lambda_step_count = 10


n_iterations = None
verbose = 1
tol = 5e-5
max_nodes=5e8

plot_rho=True
plot_A0_induced=True
save_solutions=False
read_solutions=True

inverse_lambda_step_array, lambda_critical_array = main(
        m=m,
        lambda_min=lambda_min,
        n_iterations=n_iterations,
        verbose=verbose,
        save_solutions=save_solutions,
        read_solutions=read_solutions,
        tol=tol,
        max_nodes=max_nodes,
        inverse_lambda_step_min=inverse_lambda_step_min ,
        inverse_lambda_step_max=inverse_lambda_step_max ,
        inverse_lambda_step_count=inverse_lambda_step_count ,
        )

print()
print("Save this!")
print(inverse_lambda_step_array, lambda_critical_array)


