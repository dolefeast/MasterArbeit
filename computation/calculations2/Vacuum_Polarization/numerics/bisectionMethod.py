class NoRootFoundError(Exception):
    def __init__(self, x1, x2):
        self.x1 = x1
        self.x2 = x2
        super().__init__(x1, x2)
    def __str__(self):
        return f"f({self.x1}) * f({self.x2}) > 0"

class MaxIterError(Exception):
    def __init__(self, maxiter):
        self.maxiter = maxiter
        super().__init__(maxiter)
    def __str__(self):
        return f"Number of attempts in finding a root exceeded maxiter={self.maxiter}"


def bisectionMethod(self, fun, x1, x2, tol=1e-15, maxiter=1500):
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
        raise MaxIterError(self.bisectionMaxiter)
    if fun(x1)*fun(x2) > 0:
        raise NoRootFoundError(x1, x2)

    c = (x1 + x2) / 2
    if fun(x1)*fun(c) < 0: # Means the root was at the left of c
        return self.bisectionMethod(fun, x1, c, tol=tol, maxiter=maxiter-1)
    else:
        return self.bisectionMethod(fun, c, x2, tol=tol, maxiter=maxiter-1)

def eigenvaluesUpperBound(self, previousArray):
    """
    Given an array of guesses for the roots of a function,
    returns an array of upper bound guesses
    """

    if self.maxN == 1: return [1.9 * element for element in previousArray] 
    
    positiveSolutions = [element for element in previousArray if element > 0]
    upperBound = [0] * len(positiveSolutions)

    for index, element in enumerate(positiveSolutions[:-1]):
        upperBound[index]= (element + positiveSolutions[index+1])/2

    # Some obscure math was used  to get to this value.
    upperBound[-1] = (3 * positiveSolutions[-1] - positiveSolutions[-2])/2

    return upperBound

def eigenvaluesLowerBound(self, previousArray):
    if self.maxN == 1: return [ 0.]
    
    positiveSolutions = [element for element in previousArray if element > 0]
    lowerBound = [0] * len(positiveSolutions)

    for index, element in enumerate(positiveSolutions):
        if index==0:
            continue

        lowerBound[index]= (element + positiveSolutions[index-1])/2

    return lowerBound
