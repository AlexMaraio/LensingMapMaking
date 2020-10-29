"""

Visualisations.py file which will contain the many different plots used for visualising the lensing power spectrum
across different samples all in one place

"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


class Viz:

    def __init__(self, use_latex=True):
        sns.set(font_scale=1.75, rc={'text.usetex': use_latex})

    def plot_lin_nonlin_pspec(self, linear, non_linear):
        """
        Function that plots the linear and non-linear power spectrum on the same plot

        Args:
            linear: CAMBObject class from which the data will be read in from
            non_linear: CAMBObject class from which the data will be read in from

        Returns:
            None
        """

        if linear.non_linear:
            warnings.warn('Expected the linear input object to have a linear power spectrum!')

        if not non_linear.non_linear:
            warnings.warn('Expected the non-linear input object to have a non-linear power spectrum!')

        plt.figure(figsize=(13, 7))

        plt.loglog(non_linear.ells, non_linear.get_c_ells_dict()['W1xW1'][2:], color='g', ls='-',
                   label=r'$z_1 = 0.5$')
        plt.loglog(non_linear.ells, non_linear.get_c_ells_dict()['W2xW2'][2:], color='b', ls='-',
                   label=r'$z_2 = 2$')

        plt.loglog(linear.ells, linear.get_c_ells_dict()['W1xW1'][2:], color='g', ls='--')
        plt.loglog(linear.ells, linear.get_c_ells_dict()['W2xW2'][2:], color='b', ls='--')

        plt.title('Lensing power spectrum')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        plt.legend()
        plt.tight_layout()
        plt.show(block=True)

    def plot_varying_de(self, fiducial, de):
        """
        Function that plots two lensing power spectra: the fiducial LCDM model and the varying dark-energy wCDM model

        Args:
            fiducial (CAMBObject): Class that contains the fiducial lensing power spectrum
            de (CAMBObject): Class that contains the varying dark energy power spectrum

        Returns:
            None
        """

        plt.figure(figsize=(13, 7))

        # First, plot the LCDM model
        plt.loglog(fiducial.ells, fiducial.get_c_ells_dict()['W1xW1'][2:], color='g', ls='-', label=r'$z_1;\, w=-1$')
        plt.loglog(fiducial.ells, fiducial.get_c_ells_dict()['W2xW2'][2:], color='b', ls='-', label=r'$z_2;\, w=-1$')

        # Now plot our varying dark energy model
        plt.loglog(de.ells, de.get_c_ells_dict()['W1xW1'][2:], color='g', ls='--',
                   label=r'$z_1;\, w_0 = ' + str(de.w0) + ', w_a = ' + str(de.wa) + '$')
        plt.loglog(de.ells, de.get_c_ells_dict()['W2xW2'][2:], color='b', ls='--',
                   label=r'$z_2;\, w_0 = ' + str(de.w0) + ', w_a = ' + str(de.wa) + '$')

        # Give title and axis labels
        plt.title('Lensing power spectrum')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        plt.legend(ncol=2)
        plt.tight_layout()
        plt.show(block=True)

    def plot_omega_matter_de(self, fiducial, de):
        """
        Function that plots the density ratios of matter and dark energy as a function of scale factor for two models

        Args:
            fiducial (CAMBObject): Class that contains the fiducial lensing power spectrum
            de (CAMBObject): Class that contains the varying dark energy power spectrum

        Returns:
            None
        """

        # Set up a log space for the scale factor that starts at a=1E-4 and ends at a=1 (today)
        a = np.logspace(-4, 0, 10000)

        # Obtain the density evolution with scale factor from CAMB
        rho_fid = fiducial.results.get_background_densities(a, vars=['tot', 'cdm', 'baryon', 'de'])
        rho_de = de.results.get_background_densities(a, vars=['tot', 'cdm', 'baryon', 'de'])

        fig, (ax1, ax2) = plt.subplots(2, sharex=True, figsize=(13, 10))

        # Plot the evolution of the total matter density ratio for the two models
        ax1.loglog(a, (rho_fid['cdm'] + rho_fid['baryon']) / rho_fid['tot'], color='b', ls='-', label=r'$w=-1$')
        ax1.loglog(a, (rho_de['cdm'] + rho_de['baryon']) / rho_de['tot'], color='b', ls='--',
                   label=r'$w_0 = ' + str(de.w0) + ', w_a = ' + str(de.wa) + '$')

        # Plot the evolution of the dark energy density ratio for the two models
        ax2.loglog(a, rho_fid['de'] / rho_fid['tot'], color='g', ls='-', label=r'$w=-1$')
        ax2.loglog(a, rho_de['de'] / rho_de['tot'], color='g', ls='--',
                   label=r'$w_0 = ' + str(de.w0) + ', w_a = ' + str(de.wa) + '$')

        # Add axis labels
        ax1.set_ylabel(r'$\Omega_\textrm{m}$')
        ax2.set_ylabel(r'$\Omega_\Lambda$')
        ax2.set_xlabel(r'$a$')

        # Add legends and show the plot
        ax1.legend()
        ax2.legend()
        plt.tight_layout()
        plt.show()
