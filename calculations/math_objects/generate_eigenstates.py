import numpy as np

def generate_eigenstates(self):
    # No guess eigenstate array was given
    if self.eigenstate_array is None:
        # Do we try read the solutions?
        if self.read_solutions:
            try:
                print(f"Trying to read from λ={self.lambda_value}, m={self.m}...")
                self.read_solutions_from_file()
                print(f"Found it!")
            except FileNotFoundError:
                print(
                        f"Tried reading from λ={self.lambda_value}, m={self.m} but no file was found.",
                        )
                self.read_solutions = False
                self.generate_eigenstates()
        # We don't want to read the solutions.
        # Only option left: generate them by solutions first order in λ
        else:
            print(
                        f"\n\t Generating eigenstate guesses using {self.bcs} boundary conditions"
                    )
            n = np.arange(-self.N_mode_cutoff, self.N_mode_cutoff + 1)
            self.eigenvalue_array = (
                    np.sign(n) 
                    * np.sqrt(
                        self.m ** 2
                         + (np.pi * n) ** 2
                        )
                    )
            self.eigenstate_array = [
                    self.perturbative_eigenstate(
                        self.z,
                        omega, 
                        self.m,
                        self.lambda_value
                        )
                    for omega in self.eigenvalue_array
                    ]
            self.eigenstate_gradient_array = [
                    self.perturbative_eigenstate_gradient(
                        self.z,
                        omega, 
                        self.m,
                        self.lambda_value
                        )
                    for omega in self.eigenvalue_array
                    ]
    # For completeness. Really if self.eigenstate_array was given there is not much more to be done
    elif not self.eigenstate_array is None:
        pass
