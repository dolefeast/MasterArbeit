def findRootBisection(self, fun, x1, x2, tol=1e-7, maxiter=500):
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
        return self.findRootBisection(fun, x1, c, tol=tol, maxiter=maxiter-1)
    else:
        return self.findRootBisection(fun, c, x2, tol=tol, maxiter=maxiter-1)

def bisectionMethodUpperBound(self, previousArray):
    """
    Given an array of guesses for the roots of a function,
    returns an array of upper bound guesses
    """

    if len(previousArray) == 2: return [1.9 * element for element in previousArray] 
    
    positiveSolutions = [element for element in previousArray if element > 0]
    guess = [0] * len(positiveSolutions)

    for index, element in enumerate(positiveSolutions[:-1]):
        guess[index]= (element + positiveSolutions[index+1])/2

    # Some obscure math was used  to get to this value.
    guess[-1] = (3 * positiveSolutions[-1] - positiveSolutions[-2])/2

    return [-element for element in guess[::-1]] + guess

def bisectionMethodLowerBound(self, previousArray):
    if len(previousArray) == 2: return [0., 0.]
    
    positiveSolutions = [element for element in previousArray if element > 0]
    guess = [0] * len(positiveSolutions)

    for index, element in enumerate(positiveSolutions):
        if index==0:
            continue

        guess[index]= (element + positiveSolutions[index-1])/2

    return [-element for element in guess[::-1]] + guess
