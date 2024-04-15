import regex as re 
import numpy as np

from pathlib import Path
from scripts.float_to_str import str_to_float, float_to_str

lambda_re = re.compile("lambda_\d+_\d+")
mass_re = re.compile("mass_\d+_\d+")
float_re = re.compile("\d+_\d+")

def get_lambda_value(filename):
    return str_to_float(float_re.findall(str(filename))[0])

def get_mass(filename):
    return str_to_float(float_re.findall(str(filename))[1])

def read_files(
        lambda_value=None, 
        m=None, 
        bcs='dirichlet', 
        things=[
            'eigenvalue_array',
            'eigenstate_array',
            'eigenstate_gradient_array',
            'A0_induced',
            'rho_array'
            ],
        sig_digs=3,
        ):
    """
    Given a mass, read all the files with this mass.
    Parameters:
        lambda_value=None, Read all files with this lambda value
        m=None, Read all files with this mass
        bcs='dirichlet', Read all files with thiis bcs
        things=[
            'eigenvalue',
            'eigenstate_array',
            'eigenstate_gradient_array',
            'A0_field',
            'rho'
            ], The different properties of the Vacuum_Solution class
        sig_digs=3, Significative digits to read from the files.
    Returns
        data:dict, a dictionary whose keys are the things list, and the 
                    values at each key are the Path of the file read.
                    These should be later read using np.genfromtxt.
                        TO BE CHANGED
    """
    
    
    if lambda_value is None:
        lambda_value = '*'
    else:
        lambda_value = float_to_str(lambda_value, sig_digs)

    if m is None:
        # Since there already is a * at the end of float_id
        m = ''
    else:
        m = float_to_str(m, sig_digs)
    
    file_id = f'*lambda_{lambda_value}_mass_{m}*'

    # This creates a dictionary. Its keys are the properties of the 
    # vacuum solution class. Its values are the different FILES (to be changed)
    # that contain this data. 
    data = {}
    lambda_value_array = []
    for i, thing in enumerate(things):
        data_element = sorted(
                list(
                    Path(f'saved_solutions/{bcs}/{thing}').glob(file_id)
                    ),
                key = get_lambda_value
                )
        # Since in every case the data corresponds to the same lambda value,
        # suffices to save the lambda values only for the first "column"
        if i==0:
            lambda_value_array = [
            get_lambda_value(filename) for filename in data_element
            
            ]

        data_element = [
                np.genfromtxt(filename, delimiter=",") 
                for filename in data_element
                ]

        data[thing] = data_element

    data["lambda_value"] = lambda_value_array

    return data


def max_A0(A0_perturbation_file):
    A0_perturbation = np.fromfile(
            A0_perturbation_file, 
            dtype=float,
            sep="\n"
            )
    return max(A0_perturbation)
    
