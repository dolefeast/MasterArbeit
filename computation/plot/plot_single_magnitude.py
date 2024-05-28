def plot_single_magnitude(
        m,
        a,
        lambda_value,
        data,
        ax,
        magnitude="rho",
        ):
    """
    Plots a single 'magnitude' i.e. rho, A0_induced, eigenstate_array,...
    Parameters:
        m, the mass of the magnitude of interest
        a, the interval size of the magnitude of interest
        lambda_value, the lambda_value of the magnitude of interest
        data, an dictionary with all the 'magnitudes' in every key of the dict
        ax, the axis to which to plot magnitudes
        magnitude, the magnitude in question. rho, A0_induced
    """
    import numpy as np 

    magnitudes = ['rho', 'A0_induced']
    if magnitude not in magnitudes:
        print(f'{magnitude} must be either the charge density rho or the induced electric potential A0_induced')
        return None
        
    # Reading the magnitude from the dictionary
    magnitude_data = data[magnitude] 

    try:
        np.genfromtxt(magnitude_data, delimiter=",")
    except AttributeError:
        print("The magnitude of interest is already open")

    z = np.linspace(0, 1, len(magnitude_data)) 

    ax.plot(z, magnitude_data)
