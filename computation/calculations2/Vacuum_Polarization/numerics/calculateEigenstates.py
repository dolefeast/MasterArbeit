def calculateEigenstates(self):
    if self.bcs == "dirichlet":
        # Initial values of the IVP
        initialValues = (0, 1)
        # To find the roots of phi[bcsIndex] when applying the boundary conditions
        bcsIndex = 0 
    elif self.bcs == "neumann":
        initialValues = (1, 0)
        bcsIndex = 1 
