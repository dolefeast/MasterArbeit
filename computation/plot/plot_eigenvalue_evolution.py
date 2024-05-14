def main(
        m,
        lambda_min: float=None,
        lambda_max: float=None,
        ):

    import matplotlib.pyplot as plt
    from physics import Vacuum_Polarization
    from utils.read_files_fixed_m import read_files_fixed_m

    # Reads all the data files corresponding to a mass m
    data = read_files_fixed_m(
            m,
            read_things=['eigenvalue_array'],
            )

    fig, ax = plt.subplots(
            figsize=(16, 9),
            )

    ax.plot(
            data["lambda_value"], 
            data["eigenvalue_array"],
            'b',
            )


    ax.set_ylabel("$\omega_n$")
    ax.set_xlabel("$\lambda$")

    fig.suptitle(f"Evolution of $\omega_n$ with $\lambda$, m={m}")
    plt.tight_layout()
    plt.show()
