import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font_scale=1.75, rc={'text.usetex': True})

from lib.RedshiftWindow import RedshiftWindow
from lib.CambObject import CambObject
from lib.Visualisations import Viz

import camb

if __name__ == '__main__':
    # Rough Euclid redshift window bins
    window1 = RedshiftWindow(0.2, 0.05)
    window2 = RedshiftWindow(0.49, 0.05)
    window3 = RedshiftWindow(0.62, 0.05)
    window4 = RedshiftWindow(0.74, 0.05)
    window5 = RedshiftWindow(0.85, 0.05)
    window6 = RedshiftWindow(0.96, 0.05)
    window7 = RedshiftWindow(1.08, 0.05)
    window8 = RedshiftWindow(1.24, 0.05)
    window9 = RedshiftWindow(1.45, 0.05)
    window10 = RedshiftWindow(2.04, 0.15)

    window1 = RedshiftWindow(0.5, 0.05, enable_counts=False)
    window2 = RedshiftWindow(2.0, 0.05, enable_counts=False)

    # * Initiates three models, linear LCDM, non-linear LCDM, and non-linear wCDM
    linear = CambObject('Linear', 5000, non_linear=False)
    non_linear = CambObject('Non-linear', 2000, non_linear=True)
    non_linear_de = CambObject('Non-linear-DE', 5000, non_linear=True)  # Non-linear model but now changing dark energy

    # * Sets the window functions for our three models
    linear.set_window_functions([window1, window2])
    non_linear.set_window_functions([window1, window2])
    non_linear_de.set_window_functions([window1, window2])

    non_linear.compute_c_ells()
    non_linear.plot_1x2_map_power_spectrum(key='TxT', nside=2048, use_mask=True)
    non_linear.plot_1x2_map_power_spectrum(key='W2xW2', nside=2048, use_mask=True)

    # * Determines if we want to plot the dark energy models against each other
    plot_dark_energy = True

    if plot_dark_energy:
        # Sets the dark energy parameters for our new dark energy model
        non_linear_de.set_dark_energy(w_0=-0.9, w_a=0.1)
        non_linear_de.compute_c_ells()

        # Set up a visualisation class
        viz = Viz()

        # Plot the linear and non-linear lensing power spectra
        viz.plot_lin_nonlin_pspec(linear, non_linear)

        # Plot the lensing power spectra for a varying dark energy model
        viz.plot_varying_de(non_linear, non_linear_de)

        # Plot the matter and dark energy density ratios for LCDM and wCDM models
        viz.plot_omega_matter_de(non_linear, non_linear_de)

        viz.plot_exp_gal_dens()

    # * Start of Flask running process
    non_linear.create_fields_info()
    non_linear.write_exp_gal_dist(mean_dens=45)

    non_linear.output_c_ells()
    non_linear.split_camb_output()

    # non_linear.plot_1x2_map_power_spectrum(key='W1xW1', nside=2048)

    non_linear.write_flask_config_file(n_side=2048*1)

    non_linear.set_flask_executable('~/Documents/PhD/Codes/flask/bin/flask')

    non_linear.run_flask()
    non_linear.plot_flask_output()

    # non_linear.multiple_run_flask(3000)

    non_linear.plot_multiple_run_data()
    non_linear.plot_ridge_plot()

    # !non_linear.estimate_cl_from_alm()
    # ! non_linear.plot_map_to_alm_diff()

    # ! non_linear.trim_flask_alm_output()
    # ! non_linear.use_cpp_thingy()

    # non_linear.plot_flask_output()
