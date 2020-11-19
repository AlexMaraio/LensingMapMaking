"""

Visualisations.py file which will contain the many different plots used for visualising the lensing power spectrum
across different samples all in one place

"""

import warnings
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import camb
from camb.sources import GaussianSourceWindow
import camb.model
from scipy import interpolate as interp

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
        plt.show(block=True)

        sigma_8_fid = fiducial.results.get_sigma8()
        sigma_8_de = de.results.get_sigma8()

        print('Fidicual sigma_8:', sigma_8_fid)
        print('Dark energy sigma_8:', sigma_8_de)

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

    @staticmethod
    def find_nearest_idx(array, value):
        """
        Function from StackOverflow that returns the (approximate) index for value in array

        Args:
            array (list): Input array of values
            value (float): The value that we want to find in the array

        Returns:
            The approximate index of value in array
        """
        idx = np.searchsorted(array, value, side="left")
        if idx > 0 and (idx == len(array) or math.fabs(value - array[idx - 1]) < math.fabs(value - array[idx])):
            return idx - 1  # , array[idx - 1]
        else:
            return idx  # , array[idx]

    def redshift_evo_of_scales(self):
        """
        Function to plot the redshift evolution of different scales in the lensing and matter power spectra

        Returns:
            None
        """
        # Redshift values which we will calculate the matter & lensing power spectra at
        redshifts = np.logspace(np.log10(5), np.log10(0.1), 20)

        # Dictionaries which we will store calculation results into
        matter_powers = {'z': [], 'Pk1_l': [], 'Pk1_w': [], 'Pk2_l': [], 'Pk2_w': [], 'Pk3_l': [], 'Pk3_w': []}
        lensing_powers = {'z': [], 'Cl1_l': [], 'Cl1_w': [], 'Cl2_l': [], 'Cl2_w': [], 'Cl3_l': [], 'Cl3_w': []}

        # Create two figures, one for matter power spectra and one for lensing power spectra
        fig1, ax1 = plt.subplots(figsize=(13, 7))
        fig2, ax2 = plt.subplots(figsize=(13, 7))

        # Create a matplotlib colour-map instance using the redshift values and plasma scheme
        norm = mpl.colors.LogNorm(vmin=redshifts.min(), vmax=redshifts.max())
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap='plasma')
        cmap.set_array([])

        # Here, we set ell_max to 5000.
        lmax = 5000

        # Create a CAMB class using the Lambda CDM model
        params_LCDM = camb.CAMBparams()
        params_LCDM.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112)
        params_LCDM.InitPower.set_params(As=2.1E-9, ns=0.96)
        params_LCDM.set_for_lmax(lmax, lens_potential_accuracy=1)

        params_LCDM.Want_CMB = False
        params_LCDM.NonLinear = camb.model.NonLinear_both

        # Now create a CAMB class using the wCDM model
        params_wCDM = camb.CAMBparams()
        params_wCDM.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112)
        params_wCDM.InitPower.set_params(As=2.1E-9, ns=0.96)
        params_wCDM.set_for_lmax(lmax, lens_potential_accuracy=1)

        params_wCDM.Want_CMB = False
        params_wCDM.NonLinear = camb.model.NonLinear_both

        # Set the w0 and wa values for the wCDM model
        params_wCDM.DarkEnergy = camb.dark_energy.DarkEnergyFluid(w=-0.9, wa=0.2)

        # Go through each redshift and compute power spectra at it
        for idx, z in enumerate(redshifts):
            # Set the CAMB source window and matter power spectrum variables to compute both power spectra at this z
            params_LCDM.SourceWindows = [GaussianSourceWindow(redshift=z, source_type='lensing', sigma=0.05)]
            params_LCDM.set_matter_power(redshifts=[z], kmax=10.0, silent=True)

            params_wCDM.SourceWindows = [GaussianSourceWindow(redshift=z, source_type='lensing', sigma=0.05)]
            params_wCDM.set_matter_power(redshifts=[z], kmax=10.0, silent=True)

            # Compute results for both models
            results_LCDM = camb.get_results(params_LCDM)
            results_wCDM = camb.get_results(params_wCDM)

            # Get the matter power spectrum for both the LCDM and wCDM models
            k_l, zs_l, pk_l = results_LCDM.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints=1000)
            k_w, zs_w, pk_w = results_wCDM.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints=1000)

            # Store redshift values in dictionary
            matter_powers['z'].append(z)
            lensing_powers['z'].append(z)

            # Store the LCDM matter power spectrum values for three different k values
            matter_powers['Pk1_l'].append(pk_l[0, self.find_nearest_idx(k_l, 1E-3)])
            matter_powers['Pk2_l'].append(pk_l[0, self.find_nearest_idx(k_l, 5E-2)])
            matter_powers['Pk3_l'].append(pk_l[0, self.find_nearest_idx(k_l, 2E0)])

            # Store the wCDM matter power spectrum values for three different k values
            matter_powers['Pk1_w'].append(pk_w[0, self.find_nearest_idx(k_w, 1E-3)])
            matter_powers['Pk2_w'].append(pk_w[0, self.find_nearest_idx(k_w, 5E-2)])
            matter_powers['Pk3_w'].append(pk_w[0, self.find_nearest_idx(k_w, 2E0)])

            # Now get the lensing Cl values for both models
            cls_l = results_LCDM.get_source_cls_dict()
            cls_w = results_wCDM.get_source_cls_dict()
            ls = np.arange(2, lmax + 1)

            # Store the wCDM lensing power spectrum values for three different ell values
            lensing_powers['Cl1_l'].append(cls_l['W1xW1'][20])
            lensing_powers['Cl2_l'].append(cls_l['W1xW1'][200])
            lensing_powers['Cl3_l'].append(cls_l['W1xW1'][2000])

            # Store the wCDM lensing power spectrum values for three different ell values
            lensing_powers['Cl1_w'].append(cls_w['W1xW1'][20])
            lensing_powers['Cl2_w'].append(cls_w['W1xW1'][200])
            lensing_powers['Cl3_w'].append(cls_w['W1xW1'][2000])

            # Plot the LCDM matter and lensing power spectra on their respective plots
            ax1.loglog(k_l, pk_l[0, :], c=cmap.to_rgba(z))
            ax2.loglog(ls, cls_l['W1xW1'][2:lmax + 1], c=cmap.to_rgba(z))

        # Add a colour-bar to the plots
        fig1.colorbar(cmap, label='$z$')
        fig2.colorbar(cmap, label='$z$')

        # Create colour palette for our three k/ell values
        pal = sns.cubehelix_palette(n_colors=3, rot=0.1, light=0.8, dark=0.05, gamma=1.25, hue=0.8)

        # Plot the three k values as dashed vertical lines on the matter power spectrum plot
        ax1.axvline(1E-3, ls='--', color=pal[0], lw=2)
        ax1.axvline(5E-2, ls='--', color=pal[1], lw=2)
        ax1.axvline(2E0, ls='--', color=pal[2], lw=2)

        # Plot the three ell values as dashed vertical lines on the lensing power spectrum plot
        ax2.axvline(20, ls='--', color=pal[0], lw=2)
        ax2.axvline(200, ls='--', color=pal[1], lw=2)
        ax2.axvline(2000, ls='--', color=pal[2], lw=2)

        # Give figures axis labels & titles
        ax1.set_xlabel(r'$k \,\, [h \textrm{Mpc}^{-1}]$')
        ax1.set_ylabel(r'$P(k|z) \,\, [(h^{-1} \textrm{Mpc})^3]$')
        ax1.set_title('Matter power spectrum')

        ax2.set_xlabel(r'$\ell$')
        ax2.set_ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        ax2.set_title('Lensing power spectrum')

        fig1.tight_layout()
        fig2.tight_layout()

        plt.show(block=False)

        # * Now create a plot of the matter power spectrum of the LCDM and wCDM models at our three k values
        plt.figure(figsize=(13, 7))
        plt.semilogy(matter_powers['z'], matter_powers['Pk1_l'], c=pal[0], lw=2, ls='-', label=r'$k=1 \times 10^{-3}$')
        plt.semilogy(matter_powers['z'], matter_powers['Pk2_l'], c=pal[1], lw=2, ls='-', label=r'$k=5 \times 10^{-2}$')
        plt.semilogy(matter_powers['z'], matter_powers['Pk3_l'], c=pal[2], lw=2, ls='-', label=r'$k=2 \times 10^{0}$')

        plt.semilogy(matter_powers['z'], matter_powers['Pk1_w'], c=pal[0], lw=2, ls='--')
        plt.semilogy(matter_powers['z'], matter_powers['Pk2_w'], c=pal[1], lw=2, ls='--')
        plt.semilogy(matter_powers['z'], matter_powers['Pk3_w'], c=pal[2], lw=2, ls='--')

        plt.xlabel(r'$z$')
        plt.ylabel(r'$P(k) \,\, [(h^{-1} \textrm{Mpc})^3]$')

        plt.legend()
        plt.tight_layout()

        # * Now create a plot of the lensing power spectrum of the LCDM and wCDM models at our three ell values
        plt.figure(figsize=(13, 7))
        plt.semilogy(lensing_powers['z'], lensing_powers['Cl1_l'], c=pal[0], lw=2, ls='-', label=r'$\ell = 20$')
        plt.semilogy(lensing_powers['z'], lensing_powers['Cl2_l'], c=pal[1], lw=2, ls='-', label=r'$\ell = 200$')
        plt.semilogy(lensing_powers['z'], lensing_powers['Cl3_l'], c=pal[2], lw=2, ls='-', label=r'$\ell = 2000$')

        plt.semilogy(lensing_powers['z'], lensing_powers['Cl1_w'], c=pal[0], lw=2, ls='--')
        plt.semilogy(lensing_powers['z'], lensing_powers['Cl2_w'], c=pal[1], lw=2, ls='--')
        plt.semilogy(lensing_powers['z'], lensing_powers['Cl3_w'], c=pal[2], lw=2, ls='--')

        plt.xlabel(r'$z$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')

        plt.legend()
        plt.tight_layout()

        plt.show()

    def isolate_sigma_8(self):
        """
        Function that generates the lensing and matter power spectra for a ΛCDM and wCDM model where the σ8 value
        observed at z=0 is the same between the two models.

        Returns:
            None
        """
        # The maximum l value that the lensing plots will go up to
        lmax = 5000

        # ? LambdaCDM cosmology
        params_LCDM = camb.CAMBparams()
        params_LCDM.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112)
        params_LCDM.InitPower.set_params(As=2.1E-9, ns=0.96)
        params_LCDM.set_for_lmax(lmax, lens_potential_accuracy=1)

        params_LCDM.Want_CMB = False
        params_LCDM.WantCls = False
        params_LCDM.NonLinear = camb.model.NonLinear_both

        # ? wCDM cosmology
        params_wCDM = camb.CAMBparams()
        params_wCDM.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112)
        params_wCDM.InitPower.set_params(As=2.1E-9, ns=0.96)
        params_wCDM.set_for_lmax(lmax, lens_potential_accuracy=1)

        params_wCDM.Want_CMB = False
        params_wCDM.WantCls = False
        params_wCDM.NonLinear = camb.model.NonLinear_both

        # Set the w0 and wa values for the wCDM model
        params_wCDM.DarkEnergy = camb.dark_energy.DarkEnergyFluid(w=-0.9, wa=0.2)

        # Evaluate matter power spec at z=0
        params_LCDM.set_matter_power(redshifts=[0], kmax=10.0, silent=True)
        params_wCDM.set_matter_power(redshifts=[0], kmax=10.0, silent=True)

        # List of trial As values that encompass the observed value of ~2.1E-8
        a_s_vals = np.linspace(1e-9, 3e-9, 25)

        # Lists to store the results in
        sigma_8_l_vals = []
        sigma_8_w_vals = []

        # Go through each trial As value and evaluate the matter power spectra at each to obtain σ8
        for a_s in a_s_vals:
            params_LCDM.InitPower.set_params(As=a_s, ns=0.96)
            params_wCDM.InitPower.set_params(As=a_s, ns=0.96)

            results_l = camb.get_results(params_LCDM)
            results_w = camb.get_results(params_wCDM)

            sigma_8_l = results_l.get_sigma8_0()
            sigma_8_l_vals.append(sigma_8_l)

            sigma_8_w = results_w.get_sigma8_0()
            sigma_8_w_vals.append(sigma_8_w)

        # Turn the lists into NumPy arrays
        sigma_8_l_vals = np.array(sigma_8_l_vals)
        sigma_8_w_vals = np.array(sigma_8_w_vals)

        # * Now spline the two A_s - sigma_8 relations, note that we *invert* this spline such that we can solve for
        # * the primordial A_s value given a sigma_8 value today.

        spline_l = interp.InterpolatedUnivariateSpline(sigma_8_l_vals, a_s_vals)
        spline_w = interp.InterpolatedUnivariateSpline(sigma_8_w_vals, a_s_vals)

        # The observed value of sigma_8 that we want to choose
        sigma_8_obs = 0.75

        # Plot the A_s vs sigma_8 graphs
        fig, ax = plt.subplots(figsize=(13, 7))
        ax.plot(a_s_vals, sigma_8_l_vals, lw=3, label=r'$\Lambda$CDM', color='tab:blue')
        ax.plot(a_s_vals, sigma_8_w_vals, lw=3, label='$w$CDM', color='orange')

        # Also draw on the figure the lines that correspond to our test sigma_8 value, and it's A_s values
        ax.hlines(y=sigma_8_obs, xmin=np.min(a_s_vals), xmax=spline_w(sigma_8_obs),
                  linestyles='dashed', colors='hotpink', lw=3, label=r'$\sigma_8 = 0.75$')
        ax.vlines(x=spline_l(sigma_8_obs), ymin=np.min([sigma_8_l_vals, sigma_8_w_vals]), ymax=sigma_8_obs,
                  linestyles='dashdot', colors='tab:blue', lw=3)
        ax.vlines(x=spline_w(sigma_8_obs), ymin=np.min([sigma_8_l_vals, sigma_8_w_vals]), ymax=sigma_8_obs,
                  linestyles='dashdot', colors='orange', lw=3)

        # Give title and axis labels
        ax.set_xlabel(r'$A_s$')
        ax.set_ylabel(r'$\sigma_8$')
        ax.set_title(r'Relationship between $A_s$ and $\sigma_8$ for different cosmologies')

        # Set nice axis limits
        ax.set_xlim(left=np.min(a_s_vals), right=np.max(a_s_vals))
        ax.set_ylim(bottom=np.min([sigma_8_l_vals, sigma_8_w_vals]), top=np.max([sigma_8_l_vals, sigma_8_w_vals]))

        ax.legend()
        fig.tight_layout()

        # Set the As values for our two models to ensure that they have the same sigma_8 value today
        params_LCDM.InitPower.set_params(As=spline_l(sigma_8_obs), ns=0.96)
        params_wCDM.InitPower.set_params(As=spline_w(sigma_8_obs), ns=0.96)

        # Now set the source windows for our lensing power spectrum
        params_LCDM.WantCls = True
        params_wCDM.WantCls = True
        params_LCDM.SourceWindows = [GaussianSourceWindow(redshift=0.5, source_type='lensing', sigma=0.05),
                                     GaussianSourceWindow(redshift=2, source_type='lensing', sigma=0.05)]
        params_wCDM.SourceWindows = [GaussianSourceWindow(redshift=0.5, source_type='lensing', sigma=0.05),
                                     GaussianSourceWindow(redshift=2, source_type='lensing', sigma=0.05)]

        # Now set the redshifts for the matter power spectrum
        params_LCDM.set_matter_power(redshifts=[0, 0.5, 2], kmax=10.0, silent=True)
        params_wCDM.set_matter_power(redshifts=[0, 0.5, 2], kmax=10.0, silent=True)

        # Compute results from CAMB
        results_l = camb.get_results(params_LCDM)
        results_w = camb.get_results(params_wCDM)

        # Get sigma_8 values today to ensure that they're the same between the two models
        sigma_8_l = results_l.get_sigma8_0()
        sigma_8_w = results_w.get_sigma8_0()

        print('Sigma_8 for: LCDM ', sigma_8_l, ', wCDM ', sigma_8_w)

        # Now get the lensing Cl values for both models
        cls_l = results_l.get_source_cls_dict()
        cls_w = results_w.get_source_cls_dict()
        ells = np.arange(2, lmax + 1)

        # Now get the matter power spectrum results
        k_l, zs_l, pk_l = results_l.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints=1000)
        k_w, zs_w, pk_w = results_w.get_matter_power_spectrum(minkh=1e-4, maxkh=5, npoints=1000)

        # Plot the lensing power spectrum for our two new normalised models
        fig2, ax2 = plt.subplots(figsize=(13, 7))

        ax2.loglog(ells, cls_l['W1xW1'][2:lmax + 1], c='tab:blue', ls='--', label=r'$z_1;\,\, \Lambda$CDM')
        ax2.loglog(ells, cls_w['W1xW1'][2:lmax + 1], c='orange', ls='--', label=r'$z_1;\,\, w$CDM')

        ax2.loglog(ells, cls_l['W2xW2'][2:lmax + 1], c='tab:blue', label=r'$z_2;\,\, \Lambda$CDM')
        ax2.loglog(ells, cls_w['W2xW2'][2:lmax + 1], c='orange', label=r'$z_2;\,\, w$CDM')

        ax2.set_xlabel(r'$\ell$')
        ax2.set_ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        ax2.set_title('Lensing power spectrum')

        ax2.legend(ncol=2)
        fig2.tight_layout()

        # Plot the matter power spectrum for our two new normalised models
        fig3, ax3 = plt.subplots(figsize=(13, 7))

        ax3.loglog(k_l, pk_l[1, :], c='tab:blue', ls='--', label=r'$z_1;\,\, \Lambda$CDM')
        ax3.loglog(k_w, pk_w[1, :], c='orange', ls='--', label=r'$z_1;\,\, w$CDM')

        ax3.loglog(k_l, pk_l[2, :], c='tab:blue', label=r'$z_2;\,\, \Lambda$CDM')
        ax3.loglog(k_w, pk_w[2, :], c='orange', label=r'$z_2;\,\, w$CDM')

        ax3.set_xlabel(r'$k \,\, [h \textrm{Mpc}^{-1}]$')
        ax3.set_ylabel(r'$P(k|z) \,\, [(h^{-1} \textrm{Mpc})^3]$')
        ax3.set_title('Matter power spectrum')

        ax3.legend(ncol=2)
        fig3.tight_layout()

        # Now plot the ratios of the two spectra for the lensing & matter on the same figure using two x axes
        sns.set_style('dark')
        fig4, ax4 = plt.subplots(figsize=(13, 7))

        ax4.semilogx(ells, cls_w['W1xW1'][2:lmax + 1] / cls_l['W1xW1'][2:lmax + 1],
                     c='tab:blue', ls='--', label=r'$z_1;$ lensing', lw=3)
        ax4.semilogx(ells, cls_w['W2xW2'][2:lmax + 1] / cls_l['W2xW2'][2:lmax + 1],
                     c='tab:blue', label=r'$z_2$; lensing', lw=3)

        # Create another axis that has a twin y axis for plotting k values on
        ax4_2 = ax4.twiny()

        ax4_2.semilogx(k_l, pk_w[1, :] / pk_l[1, :],
                       c='orange', ls='--', label=r'$z_1;$ matter', lw=3)
        ax4_2.semilogx(k_l, pk_w[2, :] / pk_l[2, :],
                       c='orange', label=r'$z_2;$ matter', lw=3)
        ax4_2.semilogx(k_l, pk_w[0, :] / pk_l[0, :],
                       c='darkorange', ls='-.', label=r'$z = 0;$ matter', lw=3)

        ax4.yaxis.grid(True)
        ax4.xaxis.grid(True)
        ax4.tick_params(axis='x', which='both', colors='tab:blue', labelcolor='black')
        ax4_2.tick_params(axis='x', which='both', colors='orange', labelcolor='black')

        ax4.set_xlabel(r'$\ell$')
        ax4_2.set_xlabel(r'$k \,\, [h \textrm{Mpc}^{-1}]$')
        ax4.set_ylabel(r'Ratio of $w$CDM to $\Lambda$CDM')

        fig4.legend(ncol=3, loc="upper left", bbox_to_anchor=(0, 1), bbox_transform=ax4.transAxes, handlelength=3)
        fig4.tight_layout()

        plt.show()

        # Reset the plotting style to default
        sns.set_style('darkgrid')
