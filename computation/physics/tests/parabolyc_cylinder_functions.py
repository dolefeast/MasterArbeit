import numpy as np
from scipy.special import pbdv, gamma

def calculate_coefficients(n, m, omega_n, lambda_val):
    # Define the zeros of parabolic cylinder functions corresponding to Dirichlet boundary conditions
    a_n = np.sqrt(gamma(n + 1) / (np.sqrt(np.pi) * gamma(n + 0.5)))
    b_n = -a_n
    
    return a_n, b_n

def phi_n(z, n, m, omega_n, lambda_val):
    sqrt_lambda = np.sqrt(lambda_val)
    a_n, b_n = calculate_coefficients(n, m, omega_n, lambda_val)
    
    prefactor = (1j + 1) / sqrt_lambda
    first_term = a_n * pbdv((1j * m**2) / (2 * lambda_val) - 0.5,
                                           prefactor * (omega_n + lambda_val * (z - 0.5)))
    
    prefactor = (1j - 1) / sqrt_lambda
    second_term = b_n * pbdv(- (1j * m**2) / (2 * lambda_val) - 0.5,
                                            prefactor * (omega_n + lambda_val * (z - 0.5)))
    
    return first_term + second_term

# Define parameters
n = 1  # Index of the function
m = 1  # Some other parameter
omega_n = 1  # Some other parameter
lambda_val = 1  # Some other parameter

# Define the range of z values
z_values = np.linspace(0, 1, 100)

# Calculate phi_n(z) for each z value
phi_values = phi_n(z_values, n, m, omega_n, lambda_val)

# Print or use phi_values as needed
print(phi_values)

