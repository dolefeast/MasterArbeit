def find_root_bisection(self, fun, x1, x2, tol=1e-4, maxiter=500):
    """Uses the bisection method to find the root of a function. Works recursively
    Parameters:
        fun: the callable of which to find the roots
        x1: a guess for the root
        x2: a second guess for the root
        tol=1e-4: the precision with which the roots are searched
        maxiter=500: the maximum amount of iterations before saying no solutions were found
    Returns:
        c: the solution
    """
    if abs(fun(x1)) < tol:
        return x1
    elif abs(fun(x2)) < tol:
        return x2

    if maxiter==0:
        raise RuntimeError(f"maxiter was reached")
    if fun(x1)*fun(x2) > 0:
        raise RuntimeError(f"No root was found in the interval ({x1}, {x2})")

    c = (x1 + x2) / 2
    if fun(x1)*fun(c) < 0: # Means the root was at the left of c
        return self.find_root_bisection(fun, x1, c, tol=tol, maxiter=maxiter-1)
    else:
        return self.find_root_bisection(fun, c, x2, tol=tol, maxiter=maxiter-1)

def bisection_method_upper_bound(self, previous_array):
    """
    Given an array of guesses for the roots of a function,
    returns an array of upper bound guesses
    """
    if len(previous_array) == 2: return [1.9 * element for element in previous_array] 
    
    positive_solutions = [element for element in previous_array if element > 0]
    guess = [0] * len(positive_solutions)

    for index, element in enumerate(positive_solutions[:-1]):
        guess[index]= (element + positive_solutions[index+1])/2

    # Some obscure math was used  to get to this value.
    guess[-1] = (3 * positive_solutions[-1] - positive_solutions[-2])/2

    return [-element for element in guess[::-1]] + guess

def bisection_method_lower_bound(self, previous_array):
    if len(previous_array) == 2: return [0., 0.]
    
    positive_solutions = [element for element in previous_array if element > 0]
    guess = [0] * len(positive_solutions)

    for index, element in enumerate(positive_solutions):
        if index==0:
            continue

        guess[index]= (element + positive_solutions[index-1])/2

    return [-element for element in guess[::-1]] + guess
