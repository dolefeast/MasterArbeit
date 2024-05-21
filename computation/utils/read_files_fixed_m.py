from pathlib import Path

from utils.get_lambda_value import get_lambda_value 
from utils.float_to_str import float_to_str 
from utils.open_data_array import open_data_array 

def read_files_fixed_m_a(
        m,
        a,
        bcs: str='dirichlet',
        sig_digs: int=3,
        read_things: [str]=[],
        directory="",
        ):
    """
    For a certain given mass parameter m, read all the data files.
    Parameters:
        m, the mass to read
        sig_digs: float=None, the maximum lambda_value to read from
        bcs: the boundary conditions from which to read from 
        read_things: The 'things' to read. The read things are
                    opened using np.genfromtxt. Otherwise they 
                    are left as a Path file.
        directory: inside of the saved_solutions dir, which directory to read from
    Returns: 
        data_dict, a dictionary with 
            data_dict.keys() = [
                'lambda_value', # The corresponding lambda_value for each item in the arrays
                'eigenvalue_array', # The eigenvalue array for each of the lambda_values
                'eigenstate_array', #  --      state --
                'eigenstate_gradient_array',
                'A0_induced', # Induced A0 for each of the lambda_values
                'rho', # Induced charge density for each of the lambda_values
                ]
    """
    things=[
        'eigenvalue_array',
        'eigenstate_array',
        'eigenstate_gradient_array',
        'A0_induced',
        'rho',
        ]

    m_string = float_to_str(m, sig_digs)
    a_string = float_to_str(a, sig_digs)

    file_id = f'*mass_{m_string}_a_{a_string}_E_*'

    if directory != "":
        # If directory is given, it will need an extra / before it
        directory = "/" + directory 
    
    for i, thing in enumerate(read_things):
        if thing not in things:
            print(f"Warning: The requested item '{thing}' is not a valid input. It will not be read")

    data = {}
    lambda_value_array = []
    for i, thing in enumerate(things):
        # This is suboptimal: it sorts every data array, when the ordering
        # is always the same. A more optimal approach would be to 
        # figure out what the 'sorting permutation' is, an apply it to every
        # array.
        thing_dir = f'saved_solutions{directory}/{bcs}/{thing}'

        data_element = sorted(
                list(
                    Path(thing_dir).glob(file_id)
                    ),
                key = get_lambda_value
                )

        # Defines the lambda_value_array
        # Only need to do it in only one iteration. No need
        # for more iterations
        if i==0:
            lambda_value_array = [
            get_lambda_value(filename) for filename in data_element 
            ]


        # This way it will only read what requested. 
        data[thing] = open_data_array(
                data_element,
                read_data=thing in read_things,
                )

    data["lambda_value"] = lambda_value_array

    return data

