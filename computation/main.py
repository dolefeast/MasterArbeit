from tests.increase_lambda import main

m = 0
lambda_min = .5
lambda_step = 0.25
# lambda_step = 0.01
lambda_max = 30
n_iterations = None
verbose = 0
plot = False
save_solutions = True
read_solutions = True
tol=5e-4
a=1

# for a in range(1, 11):
for a in [1, 3, 5, 8, 10, 11]:
    main(
        m=m,
        lambda_min=lambda_min,
        lambda_step=lambda_step,
        # lambda_step=0.01,
        lambda_max=lambda_max,
        n_iterations=n_iterations,
        verbose=verbose,
        plot=plot,
        save_solutions=save_solutions,
        read_solutions=read_solutions,
        tol=tol,
        a=a,
        directory="a_evolution",
    )
