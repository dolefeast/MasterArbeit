import numpy as np

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
