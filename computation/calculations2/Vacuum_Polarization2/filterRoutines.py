import numpy as np
from scipy.interpolate import CubicSpline

def extendSignal(self, x, y, paddingSize=None):

    # if paddingSize is not specified,
    # then the padding is y itself
    if paddingSize is None:
        paddingSize = len(y)
    # wrap gives periodization of the function
    yPadded = np.pad(y, paddingSize, mode="wrap")

    # need to extend x, too
    dx = x[1] - x[0]

    xPadded = np.linspace(-dx * paddingSize, 1 + dx * paddingSize, len(yPadded))
    return xPadded, yPadded

def removeNeighbourhood(
    self,
    x,
    y,
    points: (float)=(0,1),
    neighbourhoodSize: float=None,
    forceZero: bool=True,
):
    """
    Given curve (x, y) with problematic neighbourhoods around points=(x1, x2, ...), take their neighbourhood with neighbourhoodSize=neighbourhoodSize away and interpolate around it, thus smoothing the curve.
    """

    if neighbourhoodSize is None:
        # Take out the points corresponding to the convolution array
        neighbourhoodSize = (
        1 / ( self.maxN + 1 ) 
        + 2 / self.nPoints # Without this, the algorithm takes away
                            # an "open" neighbourhood (can't be open since it's discrete), 
                            # and adding it converts it into a closed one.
                            # Corresponds to adding a dz.
            )
    xList = list(x)
    yList = list(y)
    
    xArray = np.array(x)
    
    # A python array to use its addition properties
    idx = []

    # The points to take out of the array
    for p in points:
        window = np.where(abs(xArray - p) <= neighbourhoodSize)[0].tolist()
        idx += window

    idx = np.reshape(
        np.array(idx),
        -1,  # 1-d array
    )
    xList = [x for i, x in enumerate(xList) if i not in idx]
    yList = [y for i, y in enumerate(yList) if i not in idx]

    if forceZero:  # Force the removed values to go through 0
        for p in points:
            for xIndex, xValue in enumerate(xList):
                if xValue > p:
                    xList = xList[:xIndex] + [p] + xList[xIndex:]
                    yList = yList[:xIndex] + [0] + yList[xIndex:]
                    break

    return np.array(xList), np.array(yList)

def removeAndInterpolate(
        self,
    x: [float],
    y: [float],
    points=(0, 1),
    neighbourhoodSize:float=None,
    forceZero: bool=True,
):
    xRemoved, yRemoved = self.removeNeighbourhood(
        x,
        y,
        points=points,
        neighbourhoodSize=neighbourhoodSize,
        forceZero=True
    )

    interpolatedCurve = CubicSpline(
            xRemoved,
            yRemoved,
            )

    return x, interpolatedCurve(x)  # So that both are arrays

def returnTo01(self, x, y):
    idx = np.where(
        np.logical_and(
            x >= 0,
            x <= 1,
        )
    )
    return x[idx], y[idx]

def filterRho( self,rho):
    """
    Parameters:
        rho: a noisy signal (with a very particular kind of noise)
    Returns:
        rhoFiltered: the filtered signal.
    """
    # This exploits the fact that we know exactly what the 
    # noise looks like. It can be convoluted away with
    # a very specific function. Not particularly
    # relevant.

    # As long as nPoints is exactly 8 * (max_N + 1),
    # the filtering should work fine
    # As with A0, nPoints, max_N are necessarily
    # defined somewhere before the function is called.

    # Some fine tuned parameters
    if self.bcs == "dirichlet":
        edges = 2.8
        peaks = 5.7
        middle = 6
    elif self.bcs == "neumann":
        edges = 3
        peaks = 0.4
        middle = 6

    convolutingArray = [edges,  0.0, peaks,  0.0, middle,  0.0, peaks,  0.0, edges]
    # convolutingArray = [edges, 0.0, peaks, 0.0, middle, 0.0, peaks, 0.0, edges]

    # Array has to be normalized
    convolutingArray = np.array(convolutingArray) / sum(convolutingArray)

    rhoFiltered = np.convolve(rho, convolutingArray, mode="same")
 
    return rhoFiltered # This will stield yield some noise as 
                      # this is not the full smoothing algorithm.
def extendAndFilter(
        self,
    y,
    neighbourhoodSize=None,
    paddingSize=None,
    points=(0, 1),
    forceZero=True,
):
    """
    Parameters:
    x: [float], the x-array of the signal. increasing and goes from 0 to 1
    y: [float], the y-array of the signal
    filterMethod: callable, the method to be used for filtering
    neighbourhoodSize: float, the 'diameter' of the neighbourhood to be removed around points
    paddingSize=None, the size of the padding i.e. the extension of the signal. if none, paddingSize = len(y) and therefore the signal is copied both ways
    points=(0, 1), the points around which the neighbourhood is to be removed
    forceZero:bool=True, whether after removing a certain neighbourhood, we force y(points) to go through 0
    filterParameters=(, ),	 filter parameters that the filtering script may need

    1. Extends signal by paddingSize
    2. Filters the extended signal
    3. Removes the points on a neighbourhood of the boundaries (0, 1)
     3a. If forceValues=True, force the signal to go through the forced values (initially 0)
    4. Interpolates over the void area
    5. Returns the signal only in the interval (0, 1)

    """

    # First extend the signal
    if self.bcs == "dirichlet":
        x, y = self.extendSignal(self.z, y, paddingSize=paddingSize)
    elif self.bcs == "neumann":
        x = self.z

    # Second, filter it
    y = self.filterRho(y)
    y = self.filterRho(y)
    # if self.bcs == "neumann" and False:
    #     y = self.filterRho(y)

    # Third and fourth, remove the boundaries and interpolate
    if self.bcs == "dirichlet":
        x, y = self.removeAndInterpolate(
            x,
            y,
            points=points,
            neighbourhoodSize=neighbourhoodSize,
            forceZero=forceZero,
        )
        # Fifth return the signal only in the interval 0 to 1
        x, y = self.returnTo01(x, y)
    elif self.bcs == "neumann":
        x, y = self.removeNeighbourhood(
            x,
            y,
            points=points,
            neighbourhoodSize=neighbourhoodSize,
            forceZero=False,
        )

        # Linear extrapolation. Should be done in its own function.
        deriv = (y[1]-y[0])/(x[1]-x[0])
        startX = [z for z in self.z if z < x[0]]
        endX = [z for z in self.z if z > x[-1]]

        startY = [deriv * (z-x[0]) + y[0] for z in startX]
        endY = [deriv * (z-x[-1]) + y[-1] for z in endX]

        y = startY + list(y) + endY
    return x, y
