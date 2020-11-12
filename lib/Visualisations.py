"""

Visualisations.py file which will contain the many different plots used for visualising the lensing power spectrum
across different samples all in one place

"""

import warnings
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from .CambObject import CambObject


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

        if not linear.galaxy_dens == non_linear.galaxy_dens:
            raise RuntimeError('Both input CambObjects must either have galaxy counts, or neither have')

        z1 = 'W2xW2' if linear.galaxy_dens else 'W1xW1'
        z2 = 'W4xW4' if linear.galaxy_dens else 'W2xW2'

        plt.figure(figsize=(13, 7))

        plt.loglog(non_linear.ells, non_linear.get_c_ells_dict()[z1][2:], color='g', ls='-',
                   label=r'$z_1 = 0.5$')
        plt.loglog(non_linear.ells, non_linear.get_c_ells_dict()[z2][2:], color='b', ls='-',
                   label=r'$z_2 = 2$')

        plt.loglog(linear.ells, linear.get_c_ells_dict()[z1][2:], color='g', ls='--')
        plt.loglog(linear.ells, linear.get_c_ells_dict()[z2][2:], color='b', ls='--')

        plt.title('Lensing power spectrum')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        plt.legend()
        plt.tight_layout()
        plt.show(block=True)

    def plot_varying_de(self, fiducial, de):
        """
        Function that plots two power spectra: the fiducial LCDM model and the varying dark-energy wCDM model.
        Here, we plot the lensing power spectrum and the matter power spectrum.

        Args:
            fiducial (CambObject): Class that contains the fiducial lensing power spectrum
            de (CambObject): Class that contains the varying dark energy power spectrum

        Returns:
            None
        """

        if not fiducial.galaxy_dens == de.galaxy_dens:
            raise RuntimeError('Both input CambObjects must either have galaxy counts, or neither have')

        z1 = 'W2xW2' if fiducial.galaxy_dens else 'W1xW1'
        z2 = 'W4xW4' if fiducial.galaxy_dens else 'W2xW2'

        plt.figure(figsize=(13, 7))

        # First, plot the LCDM model
        plt.loglog(fiducial.ells, fiducial.get_c_ells_dict()[z1][2:], color='g', ls='-', label=r'$z_1;\, w=-1$')
        plt.loglog(fiducial.ells, fiducial.get_c_ells_dict()[z2][2:], color='b', ls='-', label=r'$z_2;\, w=-1$')

        # Now plot our varying dark energy model
        plt.loglog(de.ells, de.get_c_ells_dict()[z1][2:], color='g', ls='--',
                   label=r'$z_1;\, w_0 = ' + str(de.w0) + ', w_a = ' + str(de.wa) + '$')
        plt.loglog(de.ells, de.get_c_ells_dict()[z2][2:], color='b', ls='--',
                   label=r'$z_2;\, w_0 = ' + str(de.w0) + ', w_a = ' + str(de.wa) + '$')

        # Give title and axis labels
        plt.title('Lensing power spectrum')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        plt.legend(ncol=2)
        plt.tight_layout()
        plt.show(block=True)

        # Now, if we're doing galaxy counts, plot the comparison between the two models
        if fiducial.galaxy_dens:
            plt.figure(figsize=(13, 7))

            # First, plot the LCDM model
            plt.loglog(fiducial.ells, fiducial.get_c_ells_dict()['W1xW1'][2:], color='g', ls='-', label=r'$z_1;\, w=-1$')
            plt.loglog(fiducial.ells, fiducial.get_c_ells_dict()['W3xW3'][2:], color='b', ls='-', label=r'$z_2;\, w=-1$')

            # Now plot our varying dark energy model
            plt.loglog(de.ells, de.get_c_ells_dict()['W1xW1'][2:], color='g', ls='--',
                       label=r'$z_1;\, w_0 = ' + str(de.w0) + ', w_a = ' + str(de.wa) + '$')
            plt.loglog(de.ells, de.get_c_ells_dict()['W3xW3'][2:], color='b', ls='--',
                       label=r'$z_2;\, w_0 = ' + str(de.w0) + ', w_a = ' + str(de.wa) + '$')

            # Give title and axis labels
            plt.title('Galaxy counts power spectrum')
            plt.xlabel(r'$\ell$')
            plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
            plt.legend(ncol=2)
            plt.tight_layout()
            plt.show(block=True)

        # Get the matter power spectrum from CAMB
        k1, z1, pk1 = fiducial.results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints=10000)
        k2, z2, pk2 = de.results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints=10000)

        # Plot the two sets of matter power spectra at different redshifts
        plt.figure(figsize=(13, 7))
        plt.loglog(k1, pk1[0, :], color='g', label=r'$z=0, \, \Lambda\textrm{CDM}$')
        plt.loglog(k2, pk2[0, :], color='b', label=r'$z=0, \, w\textrm{CDM}$')

        plt.loglog(k1, pk1[1, :], color='g', ls='--', label=r'$z=0.5, \, \Lambda\textrm{CDM}$')
        plt.loglog(k2, pk2[1, :], color='b', ls='--', label=r'$z=0.5, \, w\textrm{CDM}$')

        plt.loglog(k1, pk1[2, :], color='g', ls='-.', label=r'$z=2, \, \Lambda\textrm{CDM}$')
        plt.loglog(k2, pk2[2, :], color='b', ls='-.', label=r'$z=2, \, w\textrm{CDM}$')

        plt.xlabel(r'$k \,\, [h \textrm{Mpc}^{-1}]$')
        plt.ylabel(r'$P(k|z) \,\, [(h^{-1} \textrm{Mpc})^3]$')
        plt.title('Matter power spectrum')

        plt.legend(ncol=3)
        plt.tight_layout()
        plt.show(block=False)

        # Now plot the ratio of the matter power spectra at these redshifts
        plt.figure(figsize=(13, 7))
        plt.semilogx(k1, pk2[0, :] / pk1[0, :], color='g', lw=2, label='$z=0$')
        plt.semilogx(k1, pk2[1, :] / pk1[1, :], color='purple', lw=2, label='$z=0.5$')
        plt.semilogx(k1, pk2[2, :] / pk1[2, :], color='b', lw=2, label='$z=2$')

        plt.xlabel(r'$k \,\, [h \textrm{Mpc}^{-1}]$')
        plt.ylabel(r'$P_w(k|z) / P_\Lambda(k|z)$')
        plt.title(r'Ratio of $w$CDM P(k) to $\Lambda$CDM P(k)')

        plt.legend()
        plt.tight_layout()
        plt.show(block=False)

        k1, z1, pk1 = fiducial.results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints=10000,
                                                                 var1='Weyl', var2='Weyl')
        k2, z2, pk2 = de.results.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints=10000,
                                                           var1='Weyl', var2='Weyl')

        plt.figure(figsize=(13, 7))
        plt.loglog(k1, pk1[0, :], color='g', label=r'$z=0, \, \Lambda\textrm{CDM}$')
        plt.loglog(k2, pk2[0, :], color='b', label=r'$z=0, \, w\textrm{CDM}$')

        plt.loglog(k1, pk1[1, :], color='g', ls='--', label=r'$z=0.5, \, \Lambda\textrm{CDM}$')
        plt.loglog(k2, pk2[1, :], color='b', ls='--', label=r'$z=0.5, \, w\textrm{CDM}$')

        plt.loglog(k1, pk1[2, :], color='g', ls='-.', label=r'$z=2, \, \Lambda\textrm{CDM}$')
        plt.loglog(k2, pk2[2, :], color='b', ls='-.', label=r'$z=2, \, w\textrm{CDM}$')

        plt.xlabel(r'$k \,\, [h \textrm{Mpc}^{-1}]$')
        plt.ylabel(r'$k^2 (\phi + \psi) / 2$')
        plt.title('Weyl potential')

        plt.legend(ncol=3)
        plt.tight_layout()
        plt.show(block=False)

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

    def plot_exp_gal_dens(self):
        """
        Function that plots the expected galaxy density as a function of redshift

        Returns:
            None
        """

        z_range = np.linspace(0, 2.5, 250)

        plt.figure(figsize=(13, 7))

        plt.plot(z_range, CambObject.exp_gal_dens(z_range), lw=3, label=r'$\bar{n}_g = 30; z_m=0.9$')
        plt.plot(z_range, CambObject.exp_gal_dens(z_range, mean_dens=40), lw=3, label=r'$\bar{n}_g = 40; z_m=0.9$')
        plt.plot(z_range, CambObject.exp_gal_dens(z_range, z_m=1.2), lw=3, label=r'$\bar{n}_g = 30; z_m=1.2$')

        plt.xlabel('$z$')
        plt.ylabel('Expected surface galaxy density (gal / arcmin$^2$)')

        plt.legend()
        plt.tight_layout()
        plt.show()
