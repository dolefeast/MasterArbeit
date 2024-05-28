from utils.plot_from_0_to_1 import plot_from_0_to_1
from utils.get_lambda_index import get_lambda_index

from numpy import pi

def compare_renormalization(
        a: float,
        lambda_value: float,
        data: dict,
        ax,
        ):
    """
    Compares the two renormalization techniques used. Mode sum and Hadamard. Plots them to the given axis
    Parameters:
        data: dictionary containing the rho and A0_induced quantities of interest
        lambda_value: the lambda_value for which to plot the comparison
    Returns:
        rho: The used rho
        mode_sum_rho: the rho used in ambjorn wolframs article
    """

    A0_induced_array = data['A0_induced']
    rho_array = data['rho']
    lambda_value_array = data['lambda_value']
    
    lambda_index = get_lambda_index(
            lambda_value_array,
            lambda_value,
            )

    A0_induced  = A0_induced_array[lambda_index] 
    rho  = rho_array[lambda_index] 

    z, A0_induced = plot_from_0_to_1(A0_induced)

    total_A0 = -lambda_value * (z - 1/2) + a*A0_induced

    mode_sum_rho = rho - 1/pi * total_A0

    ax.plot(z, rho, label='Hadamard renormalization')
    ax.plot(z, mode_sum_rho, label='Hadamard renormalization')

    return z, rho, mode_sum_rho
