def float_in_array(element, array, tol=1e-3):
    for num in array:
        if abs(num - element) < tol:
            return True
    return False
