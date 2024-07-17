def filter_rho(self):
    assert self.rho.ndim == 1, "self.rho.dim should be 1-dimensional array"

    if self.bcs == 'dirichlet':
        self.z, self.rho = self.total_filtering_dirichlet(self.z, self.rho)
    elif self.bcs == 'neumann':
        raise AttributeError("neumann boundary conditions have not yet been implemented")
    else:
        raise AttributeError(f"{self.bcs} are not defined")
