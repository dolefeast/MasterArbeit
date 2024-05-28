import numpy as np
def root_mean_square(y):
    return np.sqrt(np.mean(y**2))

def main(
        m=0.,
        bcs='dirichlet',
        ):

    import matplotlib.pyplot as plt
    from utils.read_files_fixed_m import read_files_fixed_m
    from utils.plot_from_0_to_1 import plot_from_0_to_1

    data = read_files_fixed_m(
            m,
            bcs,
            read_things=['eigenstate_array'],
            )

    fig, (ax_eigenstate_array, ax_norm) = plt.subplots(2)

    eigenstate_array_array = data["eigenstate_array"]
    lambda_value_array = data["lambda_value"]

    mode = len(eigenstate_array_array[0])//2 

    for i, (lambda_value, eigenstate_array) in enumerate(
            zip(
                lambda_value_array,
                eigenstate_array_array,
                )
            ):
        alpha = 0.2 + 0.8 * (lambda_value / max(lambda_value_array)) ** 2
        ax_eigenstate_array.plot(
                *plot_from_0_to_1(
                -eigenstate_array[
                    mode
                    ]
                    ),
                'b',
                alpha=alpha,
                )
        ax_norm.plot(
                lambda_value,
                root_mean_square(eigenstate_array[mode]),
                'bx',
                )

    ax_eigenstate_array.set_ylabel('$\phi_1^{(\lambda)}(z)$')
    ax_eigenstate_array.set_xlabel('$z$')

    ax_norm.set_ylabel(r'$\sqrt{\left< \phi_1^{(\lambda)}(z)^2\right>}$')
    ax_norm.set_title('$\lambda$')

    plt.show()

if __name__ == "__main__":
    main()
