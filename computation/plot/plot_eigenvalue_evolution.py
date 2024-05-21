def main(
        m,
        a,
        lambda_min: float=None,
        lambda_max: float=None,
        directory: str="",
        ):

    import matplotlib.pyplot as plt
    from physics import Vacuum_Polarization
    from utils.read_files_fixed_m import read_files_fixed_m_a

    # Reads all the data files corresponding to a mass m
    data = read_files_fixed_m_a(
            m,
            a,
            read_things=['eigenvalue_array'],
            directory=directory,
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
