def main(
        lambda_min,
        lambda_step,
        ):
    from physics import Vacuum_Polarization
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    system = Vacuum_Polarization(lambda_min)
    mode = 0

    # The guess for the KG solution
    eigenvalue = system.eigenvalue_array[mode]
    eigenstate = system.eigenstate_array[mode]
    ax.plot(eigenstate, label='Perturbative solution')

    # KG solutions without normalization
    system.calculate_eigenstates()
    eigenvalue = system.eigenvalue_array[mode]
    eigenstate = system.eigenstate_array[mode]

    ax.plot(eigenstate.sol(system.z)[0], label='Solution w/o normalizing')

    # Normalizing states
    system.normalize_eigenstates()

    # KG solutions after normalizing
    eigenvalue = system.eigenvalue_array[mode]
    eigenstate = system.eigenstate_array[mode]
    ax.plot(eigenstate, label='Normalized solution')

    print(
            'Symplectic norm of eigenstate after normalizing', 
            system.eigenstate_norm(eigenvalue, eigenstate)
            )

    ax.legend(loc='best')
    plt.show()
