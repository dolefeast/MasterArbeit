from minimal_working_example_class import Vacuum_Polarization

directory  = "ambjorn_relax_averaging"
compute = Vacuum_Polarization(max_N=1,
        lambda_min=10,
        directory=directory,
        save=False,
        save_plots=False,
        plot=True,
        show_plots=True,
        m=0,
        )

compute.full_script()
