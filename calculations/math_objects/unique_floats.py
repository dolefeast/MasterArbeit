import numpy as np

# Idea: Keeping the first float that appears may raise some
# systemic errors. Another (slower) solution would be that
# of, for every distinct float, make an array, and take
# the average of that array as the representative of
# the value in question.


def float_in_array(element, array, tol=1e-3):
    for num in array:
        if np.abs(num - element) < tol:
            return True
    return False


def unique_floats(array, precision=1e-3):
    seen = set()
    result = []

    for num in array:
        rounded_num = round(num, int(-1 * np.floor(np.log10(precision))))
        if rounded_num not in seen:
            seen.add(rounded_num)
            result.append(num)

    return result
