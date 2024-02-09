import numpy as np

def get_eigenvalues(eigenstates):
    return np.array([eigenstate.p[0] for eigenstate in eigenstates])
