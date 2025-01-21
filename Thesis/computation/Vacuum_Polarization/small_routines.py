# Some util functions
def sign(self, x):
    # Returns sign of x
    return 2 * (x>0) - 1 

def float_to_str(self, value, sig_digs=3):
    if isinstance(value, str):
        return value
    return str(
        round(
            float(
                value
            ),
            sig_digs
        )
    ).replace(".", "_")

def str_to_float(self, value:str):
    return float(
        value.replace("_", ".")
        )

def root_mean_square(self, y):
    return np.sqrt(np.mean(y**2))

def set_config_from_dict(self, config:dict):
    from scipy.interpolate import CubicSpline

    self.eigenvalues = config["eigenvalues"]
    self.eigenstates = config["eigenstates"]
    self.A0_induced = CubicSpline(self.z, config["A0_induced"])
    self.lambda_value = config["lambda_value"]
