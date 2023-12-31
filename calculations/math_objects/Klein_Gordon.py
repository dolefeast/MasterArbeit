import numpy as np

class Klein_Gordon:
    def  __init__(self):
        self.boundary_conditions ={
                'dirichlet': self.dirichlet_boundary_conditions,
                'neumann': self.neumann_boundary_conditions}

    def differential_equation(self, z, f, p):
        """The differential equation ruling the system, after assuming 
        phi_n(z, t) = exp(-i omega_n t) phi(z)
        This is the function to input in the ODE solver"""
        e = self.phi.charge
        omega_n = p[0]
        #klein_gordon = np.array((f[1], (omega_n**2 - 2*e*self.A0.value*omega_n - e**2 * self.A0.value **2)*f[0]))
        try:
            iter(self.field_strengh)
            raise ValueError("Si le quieres meter arrays al campo modifica esto idiota!")
        except TypeError:
            field = self.A0(z)
            #field = -self.field_strengh*(z -1/2) #Means I am not trying to vary 
                                                #the field (yet)

        klein_gordon = np.array((f[1], -((omega_n - e*field)**2 - self.phi.mass**2)*f[0]))
        #klein_gordon = np.array((f[1], -omega_n ** 2 * f[0]))
        return klein_gordon

    def neumann_boundary_conditions(self,ya, yb, omega_n):
        bcs = np.array((ya[1], yb[1], ya[0] - 1))
        return bcs

    def dirichlet_boundary_conditions(self, ya, yb, p):
        omega_n, = p
        bcs = np.array((ya[0], yb[0], yb[1]-1))
        return bcs

