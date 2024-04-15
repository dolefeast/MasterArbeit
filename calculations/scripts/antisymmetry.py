import numpy as np

def antisymmetry_eigenstates(array_of_arrays):
    array_of_arrays_copy = array_of_arrays.copy()
    array_of_arrays_list_copy = list(array_of_arrays_copy)
    n = len(array_of_arrays)

    for i, eigenstate in enumerate(array_of_arrays_list_copy):
        array_of_arrays_list_copy[i] = (array_of_arrays_list_copy[i] + array_of_arrays_list_copy[-1-i][::-1])/2

    return np.array(array_of_arrays_list_copy)

def antisymmetry(array_of_arrays):
    array_of_arrays_copy = array_of_arrays.copy()
    array_of_arrays_list_copy = list(array_of_arrays_copy)
    n = len(array_of_arrays)

    for i, eigenstate in enumerate(array_of_arrays_list_copy):
        array_of_arrays_copy[i] = (array_of_arrays_copy[i] - array_of_arrays[-1-i])/2

    return np.array(array_of_arrays_copy)
