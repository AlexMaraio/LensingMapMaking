from lib.RedshiftWindow import RedshiftWindow
from lib.CambObject import CambObject
from lib.Visualisations import Viz

# Import the warnings module needed to silence them
import warnings
import matplotlib as mpl

# Silence Matplotlib and Healpy warnings, which are very annoying
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=mpl.MatplotlibDeprecationWarning)

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

    window1 = RedshiftWindow(0.5, 0.05, enable_counts=True)
    window2 = RedshiftWindow(2.0, 0.05, enable_counts=True)

    # * Initiates three models, linear LCDM, non-linear LCDM, and non-linear wCDM
    linear = CambObject('Linear', 5000, 1024, non_linear=False)
    non_linear = CambObject('Non-linear', 1000, 1024, non_linear=True)
    non_linear_de = CambObject('Non-linear-DE', 5000, 1024, non_linear=True)  # Non-linear model but now changing dark energy

    # * Sets the window functions for our three models
    linear.set_window_functions([window1, window2])
    non_linear.set_window_functions([window1, window2])
    non_linear_de.set_window_functions([window1, window2])

    non_linear.compute_c_ells()

    non_linear.experimenting_with_masks()
    """
    non_linear.plot_1x2_map_power_spectrum(key='TxT', nside=2048, use_mask=True)
    non_linear.plot_1x2_map_power_spectrum(key='W2xW2', nside=2048, use_mask=True)
    """

    # * Determines if we want to plot the dark energy models against each other
    plot_dark_energy = True

    if plot_dark_energy:
        # Sets the dark energy parameters for our new dark energy model
        non_linear_de.set_dark_energy(w_0=-0.9, w_a=0.1)
        non_linear_de.compute_c_ells()

        # Set up a visualisation class
        viz = Viz()

        viz.isolate_sigma_8_for_neutrinos()
        # Plot the effects of changing the neutrino masses
        viz.effect_of_neutrino_mass()

        # Plot the linear and non-linear lensing power spectra
        viz.plot_lin_nonlin_pspec(linear, non_linear)

        # Plot the lensing power spectra for a varying dark energy model
        viz.plot_varying_de(non_linear, non_linear_de)

        # Now look at how the effects of w_0 and w_a individually affect the spectra
        viz.dark_energy_w0_w1()

        # Plot the matter and dark energy density ratios for LCDM and wCDM models
        viz.plot_omega_matter_de(non_linear, non_linear_de)

        # Plot the expected surface galaxy density distribution
        viz.plot_exp_gal_dens()

        # Plot the effects of dark energy when sigma_8 dependence is isolated
        viz.isolate_sigma_8()

        # Plot how different scales evolve in LCDM and wCDM models
        viz.redshift_evo_of_scales()

        # Plot the effects of changing the neutrino masses
        viz.effect_of_neutrino_mass()

        # Plot the effects of changing the number of massive neutrinos
        viz.neutrino_numbers()

    # * Start of Flask running process
    non_linear.create_fields_info()
    non_linear.write_exp_gal_dist(mean_dens=45)

    non_linear.output_c_ells()
    non_linear.split_camb_output()

    # non_linear.plot_1x2_map_power_spectrum(key='W1xW1', nside=2048)

    non_linear.set_multiple_masked_output('Test1')
    non_linear.write_flask_config_file()

    non_linear.set_flask_executable('~/Documents/PhD/Codes/flask/bin/flask')

    non_linear.run_flask()
    non_linear.euclid_masks()
    # ? non_linear.plot_flask_output()
    non_linear.noise_simulations()

    # non_linear.multiple_run_flask(3000)

    non_linear.plot_multiple_run_data(used_mask=True)
    non_linear.plot_ridge_plot()

    non_linear.custom_mask()
    non_linear.theta_phi_only_mask()
    non_linear.experimenting_with_masks()

    # !non_linear.estimate_cl_from_alm()
    # ! non_linear.plot_map_to_alm_diff()

    # ! non_linear.trim_flask_alm_output()
    # ! non_linear.use_cpp_thingy()

    # non_linear.plot_flask_output()
