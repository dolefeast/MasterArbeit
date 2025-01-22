import numpy as np
from minimal_working_example import float_to_str, str_to_float
import re
from pathlib import Path
from math import ceil

float_re = re.compile("\d+_\d+")

def open_data_array(
        data_array,
        read_data: bool=True,
        ):
    """
    This will be used to read all the Path objects in the read_files dictionary.
    Parameters:
        data_array, an array consisting of Path filenames.
        read_data, if False, return the original array without doing anything
    Returns:
        data_array, the read data array
    """
    if read_data:
        data_array = [
                np.genfromtxt(filename, delimiter=",")
                for filename in data_array
                ]
    return np.array(data_array)

def get_lambda_value(filename):
    # Each filename has only the mass and lambda_value parameter
    return str_to_float(float_re.findall(str(filename.name))[2])

def read_solutions_from_file(
        m,
        a,
        lambda_value,
        directory="",
        sig_digs=3,
        bcs="dirichlet",
        desired_quantity:[str]=None,
        ):
    """
    Reads the converged solutions for a certain m, a, lambda_value.
    """

    m_string = float_to_str(m, sig_digs)
    lambda_string = float_to_str(lambda_value, sig_digs)
    a_string = float_to_str(a, sig_digs)

    file_id = f"mass_{m_string}_a_{a_string}_lambda_{lambda_string}.txt"

    if directory != "": 
        directory = "/" + directory

    # Initialize the solution family dict
    solution_family = {}

    for key in desired_quantity:
        solution_family[key]  = np.genfromtxt(
            f"saved_solutions{directory}/{bcs}/{key}/{file_id}",
            dtype=float,
            delimiter="\n",
        )

    return solution_family

def downsize_unread_solution_family(
        posix_dict,
        max_lambda_density,
        ):

    if max_lambda_density is None:
        return posix_dict
    lambda_value_array = posix_dict["lambda_value"]
    lambda_density =   len(lambda_value_array) / (lambda_value_array[-1] - lambda_value_array[0])

    if lambda_density < max_lambda_density:
        # Do nothing if lambda_density is not too big,
        # or max_lambda_density was None
        return posix_dict
    
    # Now lambda_density > max_lambda_density
    density_factor = lambda_density / max_lambda_density
    # I want the len() of the arrays to be len(posix_dict) / density_factor
    new_dict = {}
    for key, value_array in posix_dict.items():
        new_dict[key] = value_array[::ceil(density_factor)]
    
    return new_dict


def get_Posix_for_quantities(
        m,
        a,
        directory="",
        sig_digs=3,
        bcs="dirichlet",
        max_lambda_density=None,
        ):

    if directory != "":
        directory = "/" + directory

    dict_keys = ['eigenvalues', 'eigenstates', 'eigenstate_gradients', 'rho', 'A0_induced']
    # Not to fill the dictionary with crap that we don't need
    
    m = float_to_str(m, sig_digs)
    a = float_to_str(a, sig_digs)

    file_regex = f"mass_{m}_a_{a}_lambda*"

    solution_family_dict = {}
    lambda_value_array = []

    for i, quantity in enumerate(dict_keys):
        data_element = sorted(
                list(
                    Path(f'saved_solutions{directory}/{bcs}/{quantity}').glob(file_regex)
                    ),
                key = get_lambda_value
                )

        solution_family_dict[quantity] = data_element

        # Enough with doing it once
        if i==0:
            lambda_value_array = [
            get_lambda_value(filename) for filename in data_element 
            ]

    solution_family_dict["lambda_value"] = lambda_value_array

    solution_family_dict = downsize_unread_solution_family(
            solution_family_dict, 
            max_lambda_density,
            )

    return solution_family_dict

def open_Posix_dict(
        posix_dict,
        desired_quantities: [str]=None,
        ):
    """
    Given a dictionary with families of solutions that store Path of solutions,
    open the desired quantities.
    Parameters:
        posix_dict: A dictionary with the keys ['eigenvalues', 'eigenstates', 'eigenstate_gradients', 'rho', 'A0_induced'],
                and items with the path of each of the files.
    """
    
    dict_keys = ['eigenvalues', 'eigenstates', 'eigenstate_gradients', 'rho', 'A0_induced']
    # Not to fill the dictionary with crap that we don't need
    final_desired_quantities = [] # To be able to ignore queries that were miswriten
    if desired_quantities is None:
        final_desired_quantities = dict_keys
    else: 
        if isinstance(desired_quantities, str):
            # Means desired quantities was a single str input
            desired_quantities = [desired_quantities]
        for quantity in desired_quantities:
            if quantity not in dict_keys:
                print(f"Solicited quantity {quantity} is not accesible. I'm skipping it.")
            else:
                final_desired_quantities.append(quantity)

    solution_family_array = {}

    for quantity in final_desired_quantities:
        data_element = open_data_array(posix_dict[quantity])
        solution_family_array[quantity] = data_element

    try:
        solution_family_array["lambda_value"] = posix_dict["lambda_value"]
    except KeyError:
        pass

    return solution_family_array


if __name__ == "__main__":
    m = 0
    a = 1
    lambda_value = 7
    directory = "ambjorn"

    posix_dict = get_Posix_for_quantities(m, a, directory=directory)
    
    solution_family_array = open_Posix_dict(posix_dict, desired_quantities=["eigenvalues"])
    print(solution_family_array)
