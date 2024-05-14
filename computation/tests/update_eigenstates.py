def main(
        m,
        lambda_value,
        ):
    from physics import Vacuum_Polarization
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(16, 9))

    system = Vacuum_Polarization(
            lambda_value=lambda_value,
            m=m,
            )

    for i in range(3):
        system.update_eigenstates()
        ax.plot(
                system.z,
                system.rho,
                alpha = (1+i)/3,
                )
    plt.show()
