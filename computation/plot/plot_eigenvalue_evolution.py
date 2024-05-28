def plot_eigenvalue_evolution(
        m,
        a,
        ax,
        directory="",
        ):

    from utils.read_files_fixed_m_a import read_files_fixed_m_a

    # Reads all the data files corresponding to a mass m
    data = read_files_fixed_m_a(
            m,
            a,
            read_things=['eigenvalue_array'],
            directory=directory,
            )

    ax.plot(
            data["lambda_value"], 
            data["eigenvalue_array"],
            'b',
            )
