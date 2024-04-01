import numpy as np
from scripts.iterable_output import iterable_output

def antisymmetry2(array):
    n = len(array)
    first_half = array[:n//2]
    second_half = array[:n//2-1:-1]

    return np.concatenate((first_half - second_half, -first_half[::-1] + second_half[::-1]))/2

def antisymmetry_eigenstates(array_of_arrays):
    array_of_arrays_copy = array_of_arrays.copy()
    array_of_arrays_list_copy = list(array_of_arrays_copy)
    n = len(array_of_arrays)

    for i, eigenstate in enumerate(array_of_arrays_list_copy):
        array_of_arrays_list_copy[i] = (array_of_arrays_list_copy[i] - array_of_arrays_list_copy[-1-i])/2

    return np.array(array_of_arrays_list_copy)

def antisymmetry(array_of_arrays):
    array_of_arrays_copy = array_of_arrays.copy()
    array_of_arrays_list_copy = list(array_of_arrays_copy)
    n = len(array_of_arrays)

    for i, eigenstate in enumerate(array_of_arrays_list_copy):
        array_of_arrays_copy[i] = (array_of_arrays_copy[i] - array_of_arrays[-1-i])/2

    return np.array(array_of_arrays_list_copy)
