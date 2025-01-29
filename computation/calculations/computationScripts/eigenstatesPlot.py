from Vacuum_Polarization import Vacuum_Polarization
import matplotlib.pyplot as plt  

def main(**kwargs):
    def fixedPointPlot(fun, axInduced, axFull):
        def inner():
            # Update
            fun()

            # Plot the new values
            axInduced[0].plot(
                    compute.z, compute.eigenstates[compute.maxN],  color=compute.color, 
                    label=f"$\lambda={compute.lambdaValue}$"
                    )

        return inner

    compute = Vacuum_Polarization(
            **kwargs
            )    

    fig1, axInduced = plt.subplots(1, 2, figsize=(12, 4))
    fig2, axFull = plt.subplots(1, 2, figsize=(12, 4))

    compute.coreIteration = fixedPointPlot(compute.coreIteration, axInduced, axFull)

    try:
        compute.fullScript()
    except KeyboardInterrupt:
        pass

    axFull[0].set_xlabel(r"$A_0^{(n)} + \frac{\lambda}{2}$")
    axFull[1].set_xlabel("n")

    axInduced[0].set_xlabel("$A_0^{ ( n ) }$ induced")
    axInduced[1].set_xlabel("n")

    axFull[0].set_ylabel(r"$A_0^{(n+1)}  + \frac{\lambda}{2}$")
    axFull[1].set_ylabel("$A_0^n$")

    axInduced[0].set_ylabel("$A_0^{(n+1)}$ induced")
    axInduced[1].set_ylabel("$A_0^n$")

    fig1.tight_layout()
    fig2.tight_layout()
    if compute.showPlots == False:
        plt.show()
