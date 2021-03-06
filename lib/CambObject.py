import os
import subprocess
import random
import pathlib
import time
import ctypes
import copy
import mpmath
import numpy as np
from scipy import stats as scistats
from scipy import constants as sciconst
from scipy import special as scispec
from scipy import interpolate as sciinterp
import pandas as pd
import functools
import camb
from camb.sources import GaussianSourceWindow
import healpy as hp
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib import colors
from matplotlib import animation as pltanim
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

from .flask_scripts.prepCambInput import split_files
from .flask_scripts.camb2info import XavierShift
from .Planck_colourmap import make_planck_colour_map


sns.set(font_scale=1.75, rc={'text.usetex': True})

num_ell = 1


class CambObject:

    def __init__(self, folder_path, lmax, n_side, non_linear=True):
        """
        General class that contains functions necessary to first compute the lensing power spectra for a given
        cosmology, save this power spectrum to the disk, input this spectra into Flask, run Flask multiple times
        to obtain multiple realisations of the Cl coefficients, and plotting the results.

        Args:
            folder_path (str): String that identifies where output data should be stored under the ./Data/ sub-directory
            lmax (int): Integer specifying the maximum l value that the lensing power spectrum should be computed to.
            n_side (int): Integer specifying the HealPix N_side parameter for any map constructed for the class.
                          Must be a power of 2 (e.g. 512, 1024, 2048, ...)
            non_linear (bool): Boolean value which indicates if we want to evaluate the non-linear
                               matter power spectrum (and thus the non-linear lensing spectrum)
        """

        # Set the path to where we will be storing data related to this class
        self.name = folder_path
        self.folder_path = 'Data/' + folder_path + '/'

        # If this folder doesn't already exist, then create it
        if not os.path.isdir(self.folder_path):
            os.makedirs(self.folder_path)

        # Store our ell_max and non_linear parameters in the class
        self.ell_max = lmax
        self.non_linear = non_linear

        # Create a vector which stores the ell values which the power spectrum will be evaluated over
        self.ells = np.arange(2, lmax + 1)

        # Are we computing galaxy number counts as well as a lensing power spectra?
        self.galaxy_dens = False

        # Create a CAMB parameters object in the class, and set the default cosmology.
        # TODO: consider changing the cosmology?
        self.params = camb.CAMBparams()
        self.params.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112)
        self.params.InitPower.set_params(As=2.1E-9, ns=0.96)
        self.params.set_for_lmax(self.ell_max, lens_potential_accuracy=1)

        # Parameters for evaluating the matter power spectra over. TODO: Think about changing redshifts?
        self.params.set_matter_power(redshifts=[0., 0.5, 2.0], kmax=10.0, silent=True)
        self.matter_power = camb.get_results(self.params)

        # Computing the CMB power spectrum is useful as it allows us to plot the TT power spectrum too
        self.params.Want_CMB = True

        # If we're doing a non-linear evaluation, then get CAMB to use non-linear models
        if self.non_linear:
            self.params.NonLinear = camb.model.NonLinear_both
            self.params.set_nonlinear_lensing(True)

            # * This can be used to change the version of HaloFit from which the non-linear spectrum is computed from
            self.params.NonLinearModel.set_params(halofit_version='mead')

        # Else, stick to the linear regime
        else:
            self.params.set_nonlinear_lensing(False)

        # Dark energy parameters
        self.w0 = -1
        self.wa = 0

        # Class members which we will store calculation results in
        self.window_functions = None
        self.num_redshift_bins = None
        self.results = None
        self.c_ells = None
        self.raw_c_ells = None
        self.camb_output_filename = None

        # Information about Flask
        self.flask_executable = None
        self.flask_executed = False

        # A list that holds a list for each redshift bin of the Cl's computed by me
        self.my_cls = []
        self.my_cls_cpp = []

        # Filepath for the mask used
        self.mask_path = None

        # Fraction of sky allowed through by the mask
        self.mask_f_sky = None

        self.masks_f_sky = []

        self.masked_cl_out_dir = None

        # Ensure that the n_side parameter is an integer power of two
        if np.log2(n_side) != int(np.log2(n_side)):
            raise RuntimeError('The n_side parameter must be an integer power of 2, e.g. 512, 2048 etc!')

        # Store the n_side parameter in the class
        self.n_side = n_side

        # File extension used when saving figures, defaults to PDF's but can be changed to PNG's if needed
        self.fig_ext = '.pdf'

    def set_dark_energy(self, w_0=-1, w_a=0):
        """
        Function to set the dark energy equation-of-state parameters. This includes the time evolution of the equation-
        of-state through the formulation
            w(a) = w0 + wa(1 − a)

        Args:
            w_0 (float): The w_0 parameter in the model
            w_a (float): The w_a parameter in the model

        Returns:
            None
        """

        # Update class parameters
        self.w0 = w_0
        self.wa = w_a

        # Update CAMB DarkEnergy model
        self.params.DarkEnergy = camb.dark_energy.DarkEnergyFluid(w=w_0, wa=w_a)

    def set_flask_executable(self, path):
        """
        Function to set the executable path for the Flask software

        Args:
            path (str): File-path that points to the Flask executable file

        Returns:
             None
        """

        # Convert the given path to an absolute path
        path = pathlib.Path(path).expanduser().resolve()

        # First check that the file exists
        if not path.is_file():
            raise RuntimeError('The given path for Flask is not correct! Please check.')

        # Then check it is an executable
        if not os.access(path, os.X_OK):
            raise RuntimeError('The path is not a valid executable file! Please check.')

        # Run Flask with no input just to ensure that it runs correctly
        command = subprocess.run(str(path), stdout=subprocess.PIPE, universal_newlines=True)
        output = command.stdout

        if 'Testing the code...' not in output:
            raise RuntimeError('The Flask test process did not run correctly! Please check.')

        # Store the Flask path in the class
        self.flask_executable = path

    def set_window_functions(self, window_functions: list) -> None:
        """
        Function that sets the class' window function properties

        Args:
            window_functions (list): A list of window functions that are specified in the RedshiftWindow.py file

        Returns:
            None
        """
        if not isinstance(window_functions, list):
            raise RuntimeError('The argument "window_functions" needs to be a list of window functions')

        self.window_functions = window_functions
        self.num_redshift_bins = len(window_functions)

        # Since we could have more than one window function type per given window function (as we can have both lensing
        # and galaxy counts), we need to manually construct a source window list and then set the class' params to this
        source_windows = []
        for window_func in window_functions:
            source_windows += window_func.construct_camb_instance()

            # If our source functions include galaxy counts, then store this in the class
            if window_func.gal_counts:
                self.galaxy_dens = True

        # Set the CAMB source windows parameters to the provided source windows
        self.params.SourceWindows = source_windows

    def set_mask(self, mask_path):
        """
        Function to set the file path for the mask that we're interested in. Also sets the sky fraction attribute
        using the provided mask

        Args:
            mask_path (str): File path for the mask that should be used

        Returns:
            None
        """

        # Convert the given path to an absolute path
        mask_path = pathlib.Path(mask_path).expanduser().resolve()

        # First check that the file exists
        if not mask_path.is_file():
            raise RuntimeError('The given path for the mask is not correct! Please check.')

        # Set the class's mask path attribute
        self.mask_path = mask_path

        # Load in the mask as a set of bool values
        print('Reading in provided mask now')
        mask = hp.read_map(mask_path, verbose=False).astype(np.bool)

        # Ensure that the mask's N_side is the same as the classes, if not then raise a warning
        if hp.get_nside(mask) != self.n_side:
            raise RuntimeWarning('The provided mask\'s N_side does not appear to match the N_side attribute for the '
                                 'class. This may or may not be a problem.')

        # Compute the fraction of sky let through by the mask
        sky_fraction = mask.sum() / mask.size

        print('Fraction of sky allowed through by the mask is {num:.2f} %'.format(num=100 * sky_fraction))

        # Fraction of sky allowed through by the mask
        self.mask_f_sky = sky_fraction

    def compute_c_ells(self):
        """
        Function that calls CAMB to evaluate the Cl values for the given model

        Returns:
            None
        """
        start_time = time.time()
        print('Running CAMB for model: ', self.name)

        self.results = camb.get_results(self.params)

        # Get both the normalised Cl values *and* raw Cl values. Store both separately
        self.c_ells = self.results.get_source_cls_dict(lmax=self.ell_max)
        self.raw_c_ells = self.results.get_cmb_unlensed_scalar_array_dict(lmax=self.ell_max, raw_cl=True,
                                                                          CMB_unit='muK')

        print('CAMB finished in {num:.2f} seconds for model {mdl}'.format(num=time.time() - start_time, mdl=self.name))

    def get_c_ells_dict(self, key=None):
        if self.c_ells is None:
            self.compute_c_ells()

        if key is None:
            return self.c_ells

        else:
            return self.c_ells[key]

    def output_c_ells(self):
        """
        Function that uses the computed Cl values by CAMB and writes the entire set of Cl values to a file. Here, the
        columns are ell, then TT etc.

        Returns:
            None
        """

        # If the Cl's are not already computed, then compute them now
        if self.raw_c_ells is None:
            self.compute_c_ells()

        filename = self.camb_output_filename = self.folder_path + 'PowerSpecCambOut.dat'

        file = open(filename, "w")
        file.write('#    L    ')
        file.write('            '.join([str(key) for key in self.raw_c_ells.keys()]))
        file.write('\n')

        for ell in self.ells:
            file.write('  ' + str(ell) + '   ')
            for key in self.raw_c_ells.keys():
                file.write(str(self.raw_c_ells[key][ell]) + '   ')
            file.write('\n')

        file.close()

    def split_camb_output(self):
        """
        Function that uses the homogenous Cl file written by "output_c_ells" and calls the external function
        "split_files" to turn these into individual files that Flask can use

        Returns:
            None
        """
        split_files(self.camb_output_filename,
                    self.camb_output_filename.split('.')[0] + '-',
                    self.folder_path + 'fields_info.ini',
                    self.num_redshift_bins,
                    self.galaxy_dens)

    def plot_1x2_map_power_spectrum(self, key, nside=2048, use_mask=False):
        """
        Function that plots the map and power spectrum of the Cl values that are given by the "key" parameter

        Args:
            key (str): String that corresponds to the key in the Cl dictionary returned by CAMB. e.g. 'TxT', 'W1xW2'.
            nside (int): N_side parameter used when constructing the maps from the Cl values
            use_mask (bool): Bool that determines if we want to apply the Planck mask to our map

        Returns:
            None
        """
        # Check if provided key is in the Cl dictionary
        if key not in self.raw_c_ells.keys():
            raise RuntimeError('The selected key "' + str(key) + '" for plotting is not in the computed array of: '
                               + str(self.raw_c_ells.keys()))

        # Obtain the Cl values using the provided key for the results dictionary. Only plot ell >= 2
        cl = self.raw_c_ells[key][2:]

        # Obtain the Planck colour map, which we will use in plots
        planck_cmap = make_planck_colour_map()

        # Use Healpy to plot the Cl's as a map
        map1 = hp.sphtfunc.synfast(cl, nside=nside, new=True, verbose=False)

        # If we're using a mask, then perform mask related actions
        if use_mask:
            # Read in the Planck mask. Note that nside must be 2048 to work with this mask

            if nside != 2048:
                raise RuntimeError('The default Planck mask has a N_side of 2048, and so the provided N_side must be '
                                   'equal to this for this to work.')

            # Read in the mask as a set of bool values
            mask = hp.read_map('./resources/existing_maps/Planck_COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits').astype(
                np.bool)

            # Compute the masked fraction, here 1 indicates no mask and 0 would be fully masked
            mask_fraction = mask.sum() / len(mask)
            print('The fraction of the mask that is allowed through is {mask_val:.3f}'.format(mask_val=mask_fraction))

            # Convert the original map into a masked map
            masked_map = hp.ma(map1)
            masked_map.mask = np.logical_not(mask)

        hp.visufunc.mollview(map1, cmap=planck_cmap, title='Map of power spectrum for ' + str(key))
        hp.graticule(verbose=False, alpha=0.8)

        # If using a mask, then plot the masked map too
        if use_mask:
            hp.visufunc.mollview(masked_map, cmap=planck_cmap,
                                 title='Map of power spectrum for ' + str(key) + 'with mask')
            hp.graticule(verbose=False, alpha=0.8)

        plt.show()

        # Now convert our map to Cl's for comparison with the original Cl's
        cl_from_map = hp.sphtfunc.anafast(map1, lmax=self.ell_max)[2:self.ell_max - 2]
        if use_mask:
            cl_with_mask = hp.sphtfunc.anafast(masked_map, lmax=self.ell_max)[2:self.ell_max - 2]
        ells = np.arange(2, self.ell_max - 2)

        # Now, we plot both sets of Cl's
        plt.figure(figsize=(13, 7))
        plt.loglog(ells, cl[:-3], label=r'$C_\ell$ input', lw=2.5)
        if use_mask:
            plt.loglog(ells, cl_from_map, label=r'$C_\ell$ from map')
        plt.title(r'Comparison of $C_{\ell}$')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$C_{\ell}$')
        plt.legend()
        plt.tight_layout()
        plt.show()

        # Plot both of the Cl's, but now with the scaling of ell(ell + 1)/2pi in place
        plt.figure(figsize=(13, 7))
        plt.loglog(ells, cl[:-3] * ells * (ells + 1) / (2 * np.pi), label=r'$C_\ell$ input', lw=2.5, color='navy')
        plt.loglog(ells, cl_from_map * ells * (ells + 1) / (2 * np.pi), label=r'$C_\ell$ from map', color='tab:blue')
        if use_mask:
            plt.loglog(ells, cl_with_mask * ells * (ells + 1) / (2 * np.pi), label=r'$C_\ell$ from map with mask',
                       color='hotpink')
            plt.loglog(ells, cl_with_mask * ells * (ells + 1) / (2 * np.pi) / mask_fraction,
                       label=r'``Corrected" $C_\ell$ from map with mask', color='tab:orange')
        plt.title(r'Comparison of $C_{\ell}$')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def create_fields_info(self):
        """
        Function that generates the fields_info.ini files that is used by Flask to determine what type and properties
        of the fields are being calculated

        Returns:
            None
        """
        filename = 'fields_info.ini'

        file = open(self.folder_path + '/' + filename, 'w')

        file.write('# File auto-generated which contains information about the redshift distributions\n\n')
        file.write('# Field number, z bin number, mean, shift, field type (1:galaxy, 2:shear), zmin, zmax\n\n')

        # Go through each provided window function
        for index, window in enumerate(self.window_functions, start=1):

            # If we have galaxy counts, then we have two fields associated with this window function: 1. the galaxy
            # counts field, and 2. the convergence field.
            if window.gal_counts:
                file.write('\t1 \t ' + str(index) + '\t0.0000\t' + '{:.4f}'.format(XavierShift(window.redshift))
                           + '\t1\t' + '{:.4f}'.format(window.redshift - window.sigma)
                           + '\t' + '{:.4f}'.format(window.redshift + window.sigma) + '\n')

                file.write('\t2 \t ' + str(index) + '\t0.0000\t' + '{:.4f}'.format(XavierShift(window.redshift))
                           + '\t2\t' + '{:.4f}'.format(window.redshift - window.sigma)
                           + '\t' + '{:.4f}'.format(window.redshift + window.sigma) + '\n')

            # Else, we just have the convergence field.
            else:
                file.write('\t1 \t ' + str(index) + '\t0.0000\t' + '{:.4f}'.format(XavierShift(window.redshift))
                           + '\t2\t' + '{:.4f}'.format(window.redshift - window.sigma)
                           + '\t' + '{:.4f}'.format(window.redshift + window.sigma) + '\n')

        file.close()

    @staticmethod
    def exp_gal_dens(z, mean_dens=30, z_m=0.9):
        """
        Function that returns the expected galaxy density at redshift z following the Euclid convention of defining
        a mean density and mean redshift

        Args:
            z (float|np.array): Redshift(s) at which to evaluate the observed galaxy density at
            mean_dens (float): The mean surface density of galaxies. Default is 30 galaxies / arcmin^2, the
                               Euclid default
            z_m (float): The mean redshift. Default is 0.9, the Euclid default

        Returns
            (float): The expected surface galaxy density
        """

        # Transform the mean redshift into z0
        z_0 = z_m / np.sqrt(2)

        # Return expected density
        return 3 / (2 * z_0) * (z / z_0) ** 2 * np.exp(-(z / z_0) ** (3 / 2))

    def write_exp_gal_dist(self, **kwargs):
        """
        Function that writes a file that contains the expected galaxy densities

        Keyword Args:
            Passes keyword arguments "mean_dens" and "z_m" into the exp_gal_dens function for varying of the
            galaxy density parameters

        Returns:
            None
        """

        # Filename that we will save the information to
        filename = 'Galaxy_dens_z_selec-f1.dat'
        filename = self.folder_path + filename

        file = open(filename, 'w')

        # Print header
        file.write('# Redshift \t Expected galaxy density distribution [gal / arcmin^2]\n')

        # A rough redshift range that encompasses all wanted redshift values
        z_range = np.linspace(0, 2.5, 100)

        # Go through each redshift value and write the redshift and galaxy density at that redshift
        for z in z_range:
            file.write('{z:.3f} \t {dens:.3f} \n'.format(z=z, dens=self.exp_gal_dens(z, **kwargs)))

        # Close file once done
        file.close()

    def write_flask_config_file(self):
        """
        Function that writes a Flask configuration file to the disk

        Returns:
            None

        """

        # Do we want Flask to run with additional inputs that use the expected galaxy number densities to add noise
        # the the lensing signal?
        with_galaxy_counts = self.galaxy_dens

        filename = 'FlaskInput.config'
        filename = self.folder_path + self.masked_cl_out_dir + '/' + filename

        file = open(filename, 'w')

        file.write('# This is a config file auto-generated by my code\n')
        file.write('# This file contains all of the cosmological, input and output data for the galaxy lensing\n\n\n')

        file.write('## Simulation basics ##\n\n')
        file.write('DIST: \t GAUSSIAN \n')
        file.write('RNDSEED: \t ' + str(random.randint(1, 10000)) + ' \n')
        file.write('POISSON: \t 1 \n')  # Enable Poisson galaxy noise

        file.write('\n## Cosmology ##\n\n')
        file.write('OMEGA_m: \t 0.3 \n')
        file.write('OMEGA_L: \t 0.7 \n')
        file.write('W_de: \t -1.0 \n\n')
        file.write('ELLIP_SIGMA: \t 0.21 \n')  # Intrinsic galaxy dispersion taken from 2010.12382
        file.write('GALDENSITY: \t 30 \n\n')  # Average surface density of galaxies for Euclid, same paper

        file.write('\n## Input data ## \n\n')
        file.write('FIELDS_INFO: \t fields_info.ini \n')
        file.write('CHOL_IN_PREFIX: \t 0 \n')
        file.write('CL_PREFIX: \t PowerSpecCambOut- \n')
        file.write('ALLOW_MISS_CL: \t 0 \n')
        file.write('SCALE_CLS: \t 1.0 \n')
        file.write('WINFUNC_SIGMA: \t -1 \n')
        file.write('APPLY_PIXWIN: \t 0 \n')  # * should this be zero?
        file.write('SUPPRESS_L: \t -1 \n')
        file.write('SUP_INDEX: \t -1 \n\n')

        file.write('\n## Survey selection functions ##\n\n')
        file.write('SELEC_SEPARABLE: \t 1 \n')
        file.write('SELEC_PREFIX: \t 0 \n')

        # Here, if we are using galaxy counts then provide Flask with the expected galaxy densities as a function
        # of redshift
        if with_galaxy_counts:
            file.write('SELEC_Z_PREFIX: \t Galaxy_dens_z_selec- \n')
        else:
            file.write('SELEC_Z_PREFIX: \t 0\n')

        file.write('SELEC_SCALE: \t 1 \n')
        file.write('SELEC_TYPE: \t 0 \n')

        # Are we using a mask?
        if self.mask_path is None:
            file.write('STARMASK: \t 0 \n\n')
        else:
            # file.write('STARMASK: \t ../../resources/Euclid_masks/Euclid-gal-mask-2048.fits \n\n')
            file.write('STARMASK: \t' + str(self.mask_path) + '\n\n')

        file.write('\n## Multipole information ##\n\n')
        file.write('EXTRAP_DIPOLE: \t 0 \n')
        file.write('LRANGE: \t 2 ' + str(self.ell_max) + '\n')
        file.write('CROP_CL: \t 0 \n')
        file.write('SHEAR_LMAX: \t' + str(self.ell_max) + '\n')
        file.write('NSIDE: \t ' + str(self.n_side) + ' \n')
        file.write('USE_HEALPIX_WGTS: \t 1 \n\n')

        file.write('\n## Covariance matrix regularisation##\n\n')
        file.write('MINDIAG_FRAC: \t 1e-12 \n')
        file.write('BADCORR_FRAC: \t 0 \n')
        file.write('REGULARIZE_METHOD: \t 1 \n')
        file.write('NEW_EVAL: \t 1e-18 \n')
        file.write('REGULARIZE_STEP: \t 0.0001 \n')
        file.write('REG_MAXSTEPS: \t 1000 \n')
        file.write('ADD_FRAC: \t 1e-10 \n')
        file.write('ZSEARCH_TOL: \t 0.0001 \n\n')

        file.write('\n## Output ##\n\n')

        # If we're using galaxy counts, then we can exit Flask once we have shear with noise, else just exit with Cl's
        if with_galaxy_counts:
            file.write('EXIT_AT: \t SHEAR_FITS_PREFIX \n')
            # file.write('EXIT_AT: \t RECOVCLS_OUT \n')
            # ! file.write('EXIT_AT: \t MAPWERFITS_PREFIX \n')
        else:
            file.write('EXIT_AT: \t RECOVCLS_OUT \n')

        file.write('FITS2TGA: \t 0 \n')
        file.write('USE_UNSEEN: \t 0 \n')
        file.write('LRANGE_OUT: \t 2 ' + str(self.ell_max) + '\n')
        file.write('MMAX_OUT: \t -1 \n')
        file.write('ANGULAR_COORD: \t 0 \n')
        file.write('DENS2KAPPA: \t 0 \n\n')

        file.write('FLIST_OUT: \t 0 \n')
        file.write('SMOOTH_CL_PREFIX: \t 0 \n')
        file.write('XIOUT_PREFIX: \t 0 \n')
        file.write('GXIOUT_PREFIX: \t 0 \n')
        file.write('GCLOUT_PREFIX: \t 0 \n')
        file.write('COVL_PREFIX: \t 0 \n')
        file.write('REG_COVL_PREFIX: \t 0 \n')
        file.write('REG_CL_PREFIX: \t 0 \n')
        file.write('CHOLESKY_PREFIX: \t 0 \n')
        file.write('AUXALM_OUT: \t 0 \n')
        file.write('RECOVAUXCLS_OUT: \t 0 \n')
        file.write('AUXMAP_OUT: \t 0 \n')
        file.write('DENS2KAPPA_STAT: \t 0 \n')

        # The MAP_OUT file is simply the field values at each coordinate on the sky after cosmic variance, but before
        # applying the mask and Poisson noise
        file.write('MAP_OUT: \t 0 \n')
        file.write('MAPFITS_PREFIX: \t 0 \n')

        # These are the recovered alm & clm values that are computed from the above map. i.e. they contain are the
        # power spectrum of the raw field values that has noise only due to cosmic variance.
        file.write('RECOVALM_OUT: \t 0 \n')
        file.write('RECOVCLS_OUT: \t Unmasked_Cls.dat \n')

        # Here, we have the shear values output the convergence and shear values at each redshift bin of the lensing
        # field, but with noise only due to cosmic variance, NOT due to the mask or Poisson galaxy sampling.
        file.write('SHEAR_ALM_PREFIX: \t 0 \n')
        file.write('SHEAR_FITS_PREFIX: \t 0 \n')  # output is (kappa, gamma1, gamma2)
        file.write('SHEAR_MAP_OUT: \t 0 \n')

        # This is the map data of each field AFTER applying cosmic variance, mask, and Poisson galaxy sampling.
        # i.e. this is the most "noisiest" output, but is the closest to the observed values and what we need to model
        if with_galaxy_counts:
            file.write('MAPWERFITS_PREFIX: \t 0 \n')
        else:
            file.write('MAPWERFITS_PREFIX: \t 0 \n')
        file.write('MAPWER_OUT: \t 0 \n')

        # Here outputs galaxy ellipticities values at each coordinate that includes galaxy noise
        file.write('ELLIPFITS_PREFIX: \t 0 \n')
        file.write('ELLIP_MAP_OUT: \t 0 \n')

        file.write('CATALOG_OUT: \t 0 \n')
        file.write('\nCATALOG_COLS: \t theta phi z kappa gamma1 gamma2 ellip1 ellip2\n')
        file.write('CAT_COL_NAMES: \t 0 \n')
        file.write('CAT32BIT: \t 0 \n')
        file.write('REDUCED_SHEAR: \t 0 \n')

        # ? Masked Cl output folder, raise an error if this has not been set
        if self.masked_cl_out_dir is None:
            # Close file first before raising error
            file.close()

            raise RuntimeError('In order to save the masked Cl data correctly, the output folder is required! Please'
                               'run the "set_multiple_masked_output" function beforehand.')

        file.write('MASKED_OUTPUT_DIR: \t' + str(self.masked_cl_out_dir) + ' \n')

        file.close()

    def run_flask(self, use_rnd=None):
        """
        Function that runs the Flask executable on an already existing Flask input .config file

        Args:
            use_rnd (bool): Flag which, when enabled, generates a new random number to be used as Flask's random number
            seed - which allows the same input .config file to generate multiple different spectra

        Returns:
            None, all data is written to the disk by Flask
        """
        # Check that the Flask executable path exists before running it
        if self.flask_executable is None:
            raise RuntimeError('Please first set the location of the Flask executable before running it!')

        if not os.path.isfile(self.folder_path + 'FlaskInput.config') and \
                not os.path.isfile(self.folder_path + self.masked_cl_out_dir + '/FlaskInput.config'):
            print('The Flask input file has not been previously generated, making one now')
            self.write_flask_config_file()

        print('\nRunning Flask')

        # Time how long it takes Flask to run
        start_time = time.time()

        # Execute Flask as a subprocess run from the shell.
        command = subprocess.run(str(self.flask_executable) + ' ' +
                                 (str(self.masked_cl_out_dir) + '/' if self.masked_cl_out_dir is not None else '') +
                                 'FlaskInput.config ' +
                                 ('RNDSEED: ' + str(random.randint(1, 1000000)) if use_rnd is not None else ''),
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                                 cwd=self.folder_path, shell=True)

        # If Flask ran successfully, save it's output as a file
        if command.returncode == 0:
            print('Flask ran successfully :)')
            print('Flask took {num:.2f} seconds to run'.format(num=time.time() - start_time))
            self.flask_executed = True

            # Write Flask's output to a file
            output_file = open(self.folder_path + 'FlaskOutput.txt', 'w')
            output_file.write(command.stdout)

            output_file.write('\nWarnings generated by Flask:\n' + command.stderr + '\n')

            output_file.close()

        # Else, there was an error somewhere, print this error and exit the program
        else:
            print('Flask did not run successfully, hopefully the error is contained below')
            print(command.stdout)
            print(command.stderr)

            raise RuntimeError('Flask did not run successfully! :(')

    def multiple_run_flask(self, num_runs, use_mask=False):
        """
        Function that runs Flask code multiple times to compute an average of the Cl values

        Args:
            num_runs (int): Number of times to run the Flask code

        Returns:
            None
        """

        # Time how long the total runs took
        start_time = time.time()

        # Lists that we will store the data-frames into after each run, one for the raw Cl values and normalised values
        raw_data_dfs = []
        data_dfs = []

        # If we are using a mask, then initiate lists that store the masked Cl values
        if use_mask:
            masked_raw_data_dfs = []
            masked_data_dfs = []

        # Run Flask the given number of times to generate a power spectra each time
        for run_num in range(num_runs):
            self.run_flask(use_rnd=True)

            # TODO: change these labels to deal with galaxy counts too!
            # Read in the new power spectrum
            cl_df = pd.read_csv(self.folder_path + 'Output-Cl.dat', sep=r'\s+')

            field_num = 2 if self.galaxy_dens else 1
            column = 'Cl-f' + str(field_num) + 'z2' + 'f' + str(field_num) + 'z2'

            # Strip any spaces out of the column names
            cl_df.columns = cl_df.columns.str.lstrip()

            # Get the ell values as a numpy array
            ells = cl_df['ell'].to_numpy()

            # Create a new data-frame that stores the raw Cl values (i.e. without ell(ell+1) normalisation)
            raw_cl_df = pd.DataFrame({'Cl': cl_df[column].to_numpy()}, index=ells)

            # Get the c_ell values normalised by ell(ell+1)/2pi, as usual
            c_ells = cl_df[column].to_numpy() * self.ells * (self.ells + 1) / (2 * np.pi)

            # Create a new dummy dataframe with out C_ell values, with an index determined by the ell values
            cl_df = pd.DataFrame({'Cl': c_ells}, index=ells)

            # Append the transpose of the current data-frames to the lists
            raw_data_dfs.append(raw_cl_df.transpose())
            data_dfs.append(cl_df.transpose())

            # If using a mask, read in the masked data too
            if use_mask:
                masked_cl_df = pd.read_csv(self.folder_path + 'Masked_Converg_Cls.dat', sep=r'\s+')

                # Strip any spaces out of the column names
                masked_cl_df.columns = masked_cl_df.columns.str.lstrip()

                # Get the ell values as a numpy array
                ells = masked_cl_df['ell'].to_numpy()

                # Normalise the masked Cl values by dividing by the sky fraction
                masked_cl_df[column] /= self.mask_f_sky

                # Create a new data-frame that stores the raw Cl values (i.e. without ell(ell+1) normalisation)
                masked_raw_cl_df = pd.DataFrame({'Cl': masked_cl_df[column].to_numpy()}, index=ells)

                # Get the c_ell values normalised by ell(ell+1)/2pi, as usual
                c_ells = masked_cl_df[column].to_numpy() * self.ells * (self.ells + 1) / (2 * np.pi)

                # Create a new dummy dataframe with out C_ell values, with an index determined by the ell values
                masked_cl_df = pd.DataFrame({'Cl': c_ells}, index=ells)

                # Append the transpose of the current data-frames to the lists
                masked_raw_data_dfs.append(masked_raw_cl_df.transpose())
                masked_data_dfs.append(masked_cl_df.transpose())

        # Concatenate the all the data frames together into a single output data frame, for both data-sets
        raw_data_df = pd.concat(raw_data_dfs, ignore_index=True)
        data_df = pd.concat(data_dfs, ignore_index=True)

        if use_mask:
            masked_raw_data_df = pd.concat(masked_raw_data_dfs, ignore_index=True)
            masked_data_df = pd.concat(masked_data_dfs, ignore_index=True)

        # Prints timing statistics
        print('The total time for the runs was {num:.3f} seconds, with an average of {numpersec:.3f} seconds/run'
              .format(num=time.time() - start_time, numpersec=(time.time() - start_time) / num_runs))

        # Save the data which can then be read in later
        raw_data_df.to_csv(self.folder_path + 'AggregateRawCls.csv', index=False)
        data_df.to_csv(self.folder_path + 'AggregateCls.csv', index=False)

        if use_mask:
            masked_raw_data_df.to_csv(self.folder_path + 'AggregateRawMaskedCls.csv', index=False)
            masked_data_df.to_csv(self.folder_path + 'AggregateMaskedCls.csv', index=False)

        # Hand off to the plotting function which plots the output data
        self.plot_multiple_run_data(used_mask=use_mask)

    def set_multiple_masked_output(self, output_folder):
        """
        Function that sets the output directory for Flask runs with multiple masks. Creates the provided folder if
        it does not exist already.

        Args:
            output_folder (str): Folder to save the temporary Cl data from Flask for each run and final aggregated
                                 Cl data

        Returns:
            None
        """

        # If the provided folder does not exist, then make it
        if not os.path.isdir(self.folder_path + output_folder):
            os.makedirs(self.folder_path + output_folder)

        # Save the folder to the class
        self.masked_cl_out_dir = output_folder

    def multiple_run_flask_with_masks(self, num_runs):
        """
        Custom function that runs Flask many times, applying ten different masks to the convergence maps to see how
        changing f_sky affects the recovered data

        Args:
            num_runs (int): The number of times that Flask should be run

        Returns:
            None
        """

        if self.masked_cl_out_dir is None:
            raise RuntimeError('In order to save the masked Cl data correctly, the output folder is required! Please'
                               'run the "set_multiple_masked_output" function beforehand.')

        # Time how long the total runs took
        start_time = time.time()

        # Dictionary of indexed by different mask numbers that we will store our the unmasked and masked data into
        data_dfs = {'Mask11': [], 'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                    'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}
        # Here, 'Mask11' refers to the unmasked values

        # Run Flask the given number of times to generate a power spectra each time
        for run_num in range(num_runs):
            self.run_flask(use_rnd=True)

            for mask_num in range(1, 12):
                cl_df = pd.read_csv(self.folder_path + self.masked_cl_out_dir + '/MaskedCls' + str(mask_num) + '.dat',
                                    sep=r'\s+')

                cl_df.columns = cl_df.columns.str.lstrip()

                # Get the ell values as a numpy array
                ells = cl_df['ell'].to_numpy()

                # Get the c_ell values normalised by ell(ell+1)/2pi, as usual
                # Also normalise by 1 / f_sky for the specific mask in question
                c_ells = cl_df['Cl-f2z2f2z2'].to_numpy() * ells * (ells + 1) / (2 * np.pi) / self.masks_f_sky[
                    mask_num - 1]

                # Create a new dummy dataframe with out C_ell values, with an index determined by the ell values
                cl_df = pd.DataFrame({'Cl': c_ells}, index=ells)

                data_dfs['Mask' + str(mask_num)].append(cl_df.transpose())

        # Save results for each mask separately
        for mask_num in range(1, 12):
            data_dfs['Mask' + str(mask_num)] = pd.concat(data_dfs['Mask' + str(mask_num)], ignore_index=True)
            data_dfs['Mask' + str(mask_num)].to_csv(self.folder_path + self.masked_cl_out_dir + '/AggregateCls_Mask' +
                                                    str(mask_num) + '.csv', index=False)

        # Prints timing statistics
        print('The total time for the runs was {num:.3f} seconds, with an average of {numpersec:.3f} seconds/run'
              .format(num=time.time() - start_time, numpersec=(time.time() - start_time) / num_runs))

        # Finally, go through each MaskedCl output file and delete them once done to tidy up
        for mask_num in range(1, 12):
            try:
                os.remove(self.folder_path + self.masked_cl_out_dir + '/MaskedCls' + str(mask_num) + '.dat')
            except OSError as err:
                print("Error: couldn't delete file {file}, error is: {error}".format(file=err.filename,
                                                                                     error=err.strerror))

    def plot_multiple_run_flask_with_masks(self):
        """
        Function to plot the data that has been obtained above for many different masks

        Returns:
            None
        """

        def moving_average(x, length):
            """
            Function that calculates the moving average of a given dataset

            Args:
                x (array): Input array of data
                length (int): The number of points to average over at each point

            Returns:
                array
            """
            return np.convolve(x, np.ones(length), 'valid') / length

        # Ensure that the folder we wish to save the figures to exists, if not then make it
        if not os.path.isdir(self.folder_path + self.masked_cl_out_dir + '/Figures'):
            os.makedirs(self.folder_path + self.masked_cl_out_dir + '/Figures')

        # Set the figure path
        fig_path = self.folder_path + self.masked_cl_out_dir + '/Figures/'

        mean_cls = {'Mask11': [], 'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                    'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        mean_cls_avg = {'Mask11': [], 'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                        'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        var_cls = {'Mask11': [], 'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                   'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        var2_cls = {'Mask11': [], 'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                    'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        skew_cls = {'Mask11': [], 'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                    'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        kurt_cls = {'Mask11': [], 'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                    'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        skew_cls_avg = {'Mask11': [], 'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                        'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        kurt_cls_avg = {'Mask11': [], 'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                        'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        # List that will hold the signal-to-noise values
        stns = []

        # Go through each mask and read in the data
        mask_max = 10
        for mask_num in range(1, mask_max):

            print(f'Mask {mask_num}', end='\t', flush=True)

            mask_key = 'Mask' + str(mask_num)

            data_df = pd.read_csv(
                self.folder_path + self.masked_cl_out_dir + '/AggregateCls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                mean_val = np.mean(c_ells)

                # Compute mean
                mean_cls[mask_key].append(mean_val)

                # Compute variance
                var_cls[mask_key].append(np.var(c_ells) / (mean_val * mean_val))
                var2_cls[mask_key].append(np.var(c_ells))

                # Compute skew and kurtosis
                skew_cls[mask_key].append(scistats.skew(np.array(c_ells), bias=True))
                kurt_cls[mask_key].append(scistats.kurtosis(np.array(c_ells), bias=True))

            # Compute the moving average of the means, skew, and kurtosis using three data points at each point
            mean_cls_avg[mask_key] = moving_average(mean_cls[mask_key] / self.c_ells['W4xW4'][2:] - 1, 3)
            skew_cls_avg[mask_key] = moving_average(skew_cls[mask_key], 3)
            kurt_cls_avg[mask_key] = moving_average(kurt_cls[mask_key], 3)

            # * Now want to estimate the signal-to-noise ratio of each map
            # To do so, we first need to evaluate the covariance matrix for each map

            # First, convert our input/"true" Cl's from CAMB and turn them into a NumPy array
            true = np.array(self.c_ells['W4xW4'][2:])

            # Range of ell values considered
            ells = np.arange(2, self.ell_max + 1)

            # The signal-to-noise is defined as the sum over all l of the Cl^2 divided by the covariance matrix at
            # (l, l), in our very basic approximation!
            stn = np.sum((np.mean(data_df).to_numpy()) ** 2 / (true ** 2 / (2 * ells + 1)))

            stns.append(stn)

        print('')

        # The figure size selected for the plots
        fig_size = [11, 6]

        # Plot the signal-to-noise as a function of f_sky
        plt.figure(figsize=fig_size)
        plt.semilogx(self.masks_f_sky, stns, 'bo')

        plt.xlabel(r'$f_\textrm{sky}$')
        plt.ylabel(r'$\left[\frac{\textsc{s}}{\textsc{n}} \right]^2$')
        plt.title('A simple signal-to-noise calculation')
        plt.tight_layout()
        plt.show()

        # norm = mpl.colors.Normalize(vmin=min(self.masks_f_sky[:9]), vmax=max(self.masks_f_sky[:9]))
        norm = mpl.colors.LogNorm(vmin=min(self.masks_f_sky[:9]), vmax=max(self.masks_f_sky[:9]))
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap='plasma')
        cmap.set_array([])

        fig1, ax1 = plt.subplots(figsize=fig_size)

        for mask_num in range(1, mask_max):
            ax1.loglog(self.ells, mean_cls['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.loglog(self.ells, self.c_ells['W4xW4'][2:], ls='--', lw=3, color='cyan', label=r'CAMB $C_\ell$')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig1.tight_layout()
        fig1.savefig(fig_path + 'PowerSpec' + self.fig_ext, bbox_inches='tight')

        # * Plot of the variance divided by the average^2 Cl
        fig2, ax2 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, mask_max):
            ax2.loglog(self.ells, var_cls['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax2.loglog(self.ells, 2 / (2 * self.ells + 1), ls='--', lw=3, color='cyan', label='Cosmic variance')

        ax2.set_xlabel(r'$\ell$')
        ax2.set_ylabel(r'$\textrm{Var}[C_\ell] / \textrm{Avg}[C_\ell]^2$')
        fig2.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig2.tight_layout()
        fig2.savefig(fig_path + 'NormVariance' + self.fig_ext, bbox_inches='tight')

        # * Plot of the variance of the Cl's divided by the Gamma-function prediction
        fig3, ax3 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, mask_max):
            ax3.loglog(self.ells, np.array(var_cls['Mask' + str(mask_num)]) / (2 / (2 * self.ells + 1)), lw=2,
                       label=r'Mask num ' + str(mask_num), c=cmap.to_rgba(self.masks_f_sky[mask_num - 1], alpha=0.8))

        ax3.set_xlabel(r'$\ell$')
        ax3.set_ylabel(r'$(\textrm{Var}[C_\ell] / \textrm{Avg}[C_\ell]^2) / \Gamma$-function prediction')
        fig3.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig3.tight_layout()
        fig3.savefig(fig_path + 'NormVarianceExp' + self.fig_ext, bbox_inches='tight')

        # * Plot of the raw variance of the Cl's
        fig3, ax3 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, mask_max):
            ax3.loglog(self.ells, var2_cls['Mask' + str(mask_num)], lw=2,
                       label=r'Mask num ' + str(mask_num), c=cmap.to_rgba(self.masks_f_sky[mask_num - 1], alpha=0.8))

        ax3.set_xlabel(r'$\ell$')
        ax3.set_ylabel(r'$\textrm{Var}[C_\ell]$')
        fig3.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig3.tight_layout()
        fig3.savefig(fig_path + 'Variance' + self.fig_ext, bbox_inches='tight')

        # * Plot of the relative difference between the average Cl and input values
        fig3, ax3 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, mask_max):
            ax3.semilogx(self.ells, np.array(mean_cls['Mask' + str(mask_num)]) / self.c_ells['W4xW4'][2:] - 1,
                         lw=2, label=r'Mask num ' + str(mask_num),
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax3.set_xlabel(r'$\ell$')
        ax3.set_ylabel(r'$C_\ell / C_\ell^\textrm{In} - 1$')
        plt.yscale('symlog', linthreshy=0.0075)
        ax3.set_title(r'Relative difference between average recovered $C_\ell$ and input values')
        fig3.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig3.tight_layout()
        fig3.savefig(fig_path + 'ClRelDiff' + self.fig_ext, bbox_inches='tight')

        # * Plot of the moving average of the relative difference
        fig3, ax3 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, mask_max):
            ax3.semilogx(self.ells[1:-1], mean_cls_avg['Mask' + str(mask_num)],
                         lw=2, label=r'Mask num ' + str(mask_num),
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax3.set_xlabel(r'$\ell$')
        ax3.set_ylabel(r'Rolling average of $C_\ell / C_\ell^\textrm{In} - 1$')
        plt.yscale('symlog', linthreshy=0.0075)
        ax3.set_title(r'Relative difference between average recovered $C_\ell$ and input values')
        fig3.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig3.tight_layout()
        fig3.savefig(fig_path + 'ClRelDiff_MovAvg' + self.fig_ext, bbox_inches='tight')

        # * Plot of the sum of the differences for different f_sky
        fig3, ax3 = plt.subplots(figsize=fig_size)
        for idx, mask_num in enumerate(range(1, mask_max)):
            ax3.plot(self.masks_f_sky[idx],
                     np.sum(np.array(mean_cls['Mask' + str(mask_num)]) / self.c_ells['W4xW4'][2:] - 1), 'x',
                     lw=2, label=r'Mask num ' + str(mask_num),
                     c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax3.set_xlabel(r'$f_\textrm{sky}$')
        ax3.set_ylabel(r'$\Sigma_\ell \left[ \,\, \hat{C}_\ell / C_\ell - 1 \right]$')
        ax3.set_title(r'Sum of the relative differences in $C_\ell$ values for different $f_\textrm{sky}$')
        fig3.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig3.tight_layout()
        fig3.savefig(fig_path + 'SumRelDiff' + self.fig_ext, bbox_inches='tight')

        # * Plot of the sum of the squared-residuals for different f_sky
        fig3, ax3 = plt.subplots(figsize=fig_size)
        for idx, mask_num in enumerate(range(1, mask_max)):
            ax3.plot(self.masks_f_sky[idx],
                     np.sum((np.array(mean_cls['Mask' + str(mask_num)]) / self.c_ells['W4xW4'][2:] - 1) ** 2), 'x',
                     lw=2, label=r'Mask num ' + str(mask_num),
                     c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax3.set_xlabel(r'$f_\textrm{sky}$')
        ax3.set_ylabel(r'$\Sigma_\ell \left[ \,\, \hat{C}_\ell / C_\ell - 1 \right]^2$')
        ax3.set_title(r'Sum of the relative differences in $C_\ell$ values for different $f_\textrm{sky}$')
        fig3.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig3.tight_layout()
        fig3.savefig(fig_path + 'SumRelDiffSq' + self.fig_ext, bbox_inches='tight')

        # * Plot of the skew
        fig4, ax4 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, mask_max):
            ax4.semilogx(self.ells, skew_cls['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        k = (2 * self.ells + 1) / 2
        ax4.semilogx(self.ells, 2 / np.sqrt(k), color='cyan', ls='--', lw=3, label=r'$\Gamma$ function prediction')

        ax4.set_xlabel(r'$\ell$')
        ax4.set_ylabel(r'Skew')
        fig4.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig4.savefig(fig_path + 'Skew' + self.fig_ext, bbox_inches='tight')

        # * Plot of the rolling average of the skew
        fig4, ax4 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, mask_max):
            ax4.semilogx(self.ells[1:-1], skew_cls_avg['Mask' + str(mask_num)], lw=2,
                         label=r'Mask num ' + str(mask_num),
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        k = (2 * self.ells + 1) / 2
        ax4.semilogx(self.ells, 2 / np.sqrt(k), color='cyan', ls='--', lw=3, label=r'$\Gamma$ function prediction')

        ax4.set_xlabel(r'$\ell$')
        ax4.set_ylabel(r'Rolling average of skew')
        fig4.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig4.tight_layout()
        fig4.savefig(fig_path + 'Skew_MovAvg' + self.fig_ext, bbox_inches='tight')

        # * Plot of the kurtosis
        fig5, ax5 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, mask_max):
            ax5.semilogx(self.ells, kurt_cls['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        plt.semilogx(self.ells, 6 / k, color='cyan', ls='--', lw=3, label=r'$\Gamma$ function prediction')

        ax5.set_xlabel(r'$\ell$')
        ax5.set_ylabel(r'Excess kurtosis')
        fig5.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig5.tight_layout()
        fig5.savefig(fig_path + 'Kurt' + self.fig_ext, bbox_inches='tight')

        # * Plot of the moving average of the kurtosis
        fig5, ax5 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, mask_max):
            ax5.semilogx(self.ells[1:-1], kurt_cls_avg['Mask' + str(mask_num)], lw=2,
                         label=r'Mask num ' + str(mask_num),
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        plt.semilogx(self.ells, 6 / k, color='cyan', ls='--', lw=3, label=r'$\Gamma$ function prediction')

        ax5.set_xlabel(r'$\ell$')
        ax5.set_ylabel(r'Rolling averaged excess kurtosis')
        fig5.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig5.tight_layout()
        fig5.savefig(fig_path + 'Kurt_MovAvg' + self.fig_ext, bbox_inches='tight')

        plt.show()

    def plot_correlation_animation(self, ell_max=125):
        """
        Function that reads in previously saved Cl data for different masks and produces an animation showing how the
        correlation matrix changes with masks.

        Args:
            ell_max (int): Maximum ell value that the covariance matrix should be computed up to. 125 is a good default,
                           max value is probably around 250 or so.

        Returns:
            None
        """

        # Here, the data variable is a list of DataFrames for different masks.
        data = []

        # Go through each mask
        for idx in range(1, 10):
            print('Mask ' + str(idx), end='\t', flush=True)

            # Read in the set of Cl data
            data_df = pd.read_csv(self.folder_path + self.masked_cl_out_dir + '/AggregateCls_Mask' + str(idx) + '.csv')

            # Set up a dictionary that will have ell1, ell2, and correlation coefficient data
            pearson_data = {'x': [], 'y': [], 'r': []}

            # Go through each Cl combination and work out the Pearson correlation coefficient
            for ell1, c_ells1 in data_df.items():
                # Only go up to ell1 of ell_max which was provided
                if int(ell1) > ell_max:
                    continue

                # Only want to compute samples where ell2 < ell1
                for ell2, c_ells2 in data_df.items():
                    if int(ell2) > ell_max:  # >= int(ell1):
                        continue

                    # Use the Pearson test to get correlation coefficient
                    r_val, p_val = scistats.pearsonr(c_ells1, c_ells2)

                    # Store the calculation results in the Pearson data dictionary
                    pearson_data['x'].append(int(ell1))
                    pearson_data['y'].append(int(ell2))
                    pearson_data['r'].append(r_val)

            # Turn our dictionary into a data-frame
            pearson_df = pd.DataFrame(pearson_data)

            # Pivot the data-frame to get it in the right format for plotting
            pearson_df = pearson_df.pivot('x', 'y', 'r')

            # Append the current DataFrame to our lists
            data.append(pearson_df)

        print('')

        # For the colour-map, find the minimum and maximum correlation
        min_val = np.min([np.min(data_iter) for data_iter in data])
        max_val = np.max([np.max(data_iter) for data_iter in data])

        # Print this data
        print('Most negative correlation: {corr:.3f}'.format(corr=min_val))
        print('Most positive correlation: {corr:.3f}'.format(corr=max_val))

        # Set up the figure to start with: a heatmap with colour-bar and axis labels
        fig = plt.figure()
        sns.heatmap(data[0], vmax=1, vmin=min_val, square=True, cmap='seismic', center=0,
                    norm=colors.SymLogNorm(linthresh=0.1, vmin=min_val, vmax=1),
                    xticklabels='l2', yticklabels='l1')

        def animate(frame_num):
            """
                Function that gets called every time a new frame is drawn in the animation

            Args:
                frame_num (int): The current frame number, provided by matplotlib's FuncAnimation function. Starts at 0

            Returns:
                None
            """
            # Update the heatmap by changing the data
            sns.heatmap(data[frame_num], vmax=1, vmin=min_val, square=True, cbar=False, cmap='seismic', center=0,
                        norm=colors.SymLogNorm(linthresh=0.1, vmin=min_val, vmax=1))

            # Also update the plot title with mask number and f_sky
            plt.title('Mask ' + str(frame_num + 1) + r', fsky: {fsky:.2f} \%'.format(
                fsky=self.masks_f_sky[frame_num] * 100))

            plt.xlabel('l1')
            plt.ylabel('l2')

        # Use the maptlotlib animation module to create the animation
        anim = pltanim.FuncAnimation(fig, animate, frames=len(data), interval=750, repeat=True)

        # Create a writer class to save the animation as a GIF
        gif_writer = pltanim.PillowWriter(fps=3)

        # Now save the animation
        anim.save(self.folder_path + self.masked_cl_out_dir + '/CorrelationAnim.gif', writer=gif_writer)

        # Also show the animation to the screen
        plt.show()

    def plot_covariance_matrix(self):
        """
        Function to plot the theoretical and numerical covariance matrices for a set of masks to compare how accurate
        the Pseudo-Cl estimate of the covariance matrix is

        Returns:
            None
        """

        # The l_max that that we want to calculate the covariance matrices up to
        lmax = 1000
        ells = np.arange(2, lmax + 1)

        # Lists to store data in
        mask_cls = []
        covs = []
        covs_th = []
        num_vars = []

        # Go through each mask
        for idx in range(1, 11):
            print(f'Mask {idx}', end='\t', flush=True)

            if idx < 10:
                mask = hp.read_map(self.folder_path + f'Masks/Mask{idx}.fits', verbose=False)
            else:
                # Mask 10 is the Euclid mask
                mask = hp.read_map('./resources/Euclid_masks/Euclid-gal-mask-1024.fits', verbose=False).astype(np.bool)

            # Compute the Cl values for the mask, up to 2 * lmax
            mask_cl = hp.anafast(mask, lmax=2 * lmax)
            mask_cls.append(mask_cl[2:])

            # Use CAMB to compute the mixing matrix
            mask_matrix = camb.mathutils.scalar_coupling_matrix(mask_cl, lmax=lmax)[2:, 2:]

            # Turn the mixing matrix into the symmetric form
            mask_matrix /= (2 * ells + 1)[np.newaxis, :]

            cov_th = 2 * self.c_ells['W4xW4'][np.newaxis, 2:lmax + 1] * self.c_ells['W4xW4'][2:lmax + 1, np.newaxis] \
                     * mask_matrix

            covs_th.append(cov_th)

            var_cls = []

            data_df = pd.read_csv(self.folder_path + self.masked_cl_out_dir + '/AggregateCls_Mask' + str(idx) + '.csv')

            data_df *= self.masks_f_sky[idx - 1]

            for label, c_ells in data_df.items():
                var_cls.append(np.var(c_ells))

            num_vars.append(var_cls)

            cov = np.cov(data_df.to_numpy(), rowvar=False)

            covs.append(cov)

        # Turn our lists into an array
        num_vars = np.array(num_vars)
        print('')

        # * First, plot the Cl's for the various masks & Euclid mask

        # Read in the Euclid mask
        euclid_mask = hp.read_map('./resources/Euclid_masks/Euclid-gal-mask-1024.fits', verbose=False).astype(
            np.bool)

        # Calculate the Cl values for the mask
        euclid_mask_cls = hp.anafast(euclid_mask, lmax=2 * lmax)[2:]

        # Plot the Cl values for our masks and the Euclid mask
        plt.figure(figsize=(11, 6))
        cl_ells = np.arange(2, 2 * lmax + 1)

        # Some nice colours to plot with
        colours = ['tab:blue', 'orange', 'hotpink', 'purple', 'tab:cyan']

        # Only plot half of the masks for plot clarity
        for idx, cls in enumerate(mask_cls):
            if idx % 2 != 0:
                continue

            # Plot only even ell
            plt.loglog(cl_ells[::2], cl_ells[::2] * (cl_ells[::2] + 1) * cls[::2] / (2 * np.pi),
                       lw=1.5, c=colours[idx // 2],
                       label=r'$f_\textrm{{sky}} = {fsky:.2f} \%$'.format(fsky=self.masks_f_sky[idx] * 100))

        plt.loglog(cl_ells, cl_ells * (cl_ells + 1) * euclid_mask_cls / (2 * np.pi),
                   lw=1.5, c='lime', label='Euclid mask')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_{\ell} / 2 \pi$')
        plt.title(r'Power spectrum of mask for even $\ell$ only')
        plt.legend()
        plt.tight_layout()
        plt.show()

        # * Now plot the numerical and theoretical form of the covariance matrices

        # Here, use the MASTER code to compute it too
        import pymaster as mst

        # N_side parameter for the map to generate
        nside = 1024

        # Compute f_sky for Euclid mask
        f_sky = euclid_mask.sum() / euclid_mask.size

        # Obtain our TT convergence Cl values
        cl_tt = self.raw_c_ells['W4xW4']

        # Generate a map from our TT cl's
        map_t = hp.synfast(cl_tt, nside, verbose=False)

        # Mask the map
        map_t = hp.ma(map_t)
        map_t.mask = np.logical_not(euclid_mask)

        # The l max that we want to recover the masked map up to (MASTER requires lmax >= Nside)
        lmax = 3 * nside + 1
        cl_tt = hp.anafast(map_t, lmax=lmax)

        # The ell range for our Cl values
        ells = np.arange(0, lmax + 1)

        # Normalise the recovered Cl values
        cl_tt = ells * (ells + 1) * cl_tt / (2 * np.pi)

        # Implement 1/f_sky factor
        cl_tt /= f_sky

        # Use the MASTER code to compute the fields
        field = mst.NmtField(euclid_mask, [map_t])

        # We only use one l-mode per bin for the moment, to compare with previous results
        bins = mst.NmtBin.from_nside_linear(nside, 1)

        # Create a workspace, and compute the coupling matrix for the mask
        work = mst.NmtWorkspace()
        work.compute_coupling_matrix(field, field, bins)

        cowork = mst.NmtCovarianceWorkspace()

        cowork.compute_coupling_coefficients(field, field, field, field)

        # Compute the Gaussian covariance matrix for our mask
        covar_00_00 = mst.gaussian_covariance(cowork,
                                              0, 0, 0, 0,  # Spins of the 4 fields
                                              [cl_tt], [cl_tt], [cl_tt], [cl_tt],
                                              work, wb=work)

        # Re-define lmax here, to compare with previous results
        lmax = 1000

        # First plot just the MASTER covariance matrix
        plt.imshow(covar_00_00[0:lmax - 1, 0:lmax - 1], cmap='inferno')
        plt.xlabel(r'$\ell_1$')
        plt.ylabel(r'$\ell_2$')
        plt.title('MASTER covariance matrix')
        plt.colorbar()

        # Plot absolute difference
        fig3, ax3 = plt.subplots()
        im3 = ax3.imshow(covar_00_00[0:lmax - 1, 0:lmax - 1] - covs[9], cmap=make_planck_colour_map())
        ax3.set_xlabel(r'$\ell_1$')
        ax3.set_ylabel(r'$\ell_2$')
        ax3.set_title(f'Difference between MASTER and numerical covariance')
        fig3.colorbar(im3, label='MASTER - numerical')

        # Plot relative difference on a log scale
        fig3, ax3 = plt.subplots()
        im3 = ax3.imshow(np.log10(np.abs(covar_00_00[0:lmax - 1, 0:lmax - 1] / covs[9])), cmap=make_planck_colour_map())
        ax3.set_xlabel(r'$\ell_1$')
        ax3.set_ylabel(r'$\ell_2$')
        ax3.set_title(f'Relative difference between MASTER and numerical covariance matrices')
        fig3.colorbar(im3, label=r'$\log_{10}(\textrm{Abs}[$MASTER / numerical$])$')

        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(13, 10), sharey=True)

        # The selected mask number (ranges from 0 to 9, with 9 being Euclid)
        mask_idx = 9

        min_val = np.min([covs[mask_idx], covs_th[mask_idx]])
        max_val = np.max([covs[mask_idx], covs_th[mask_idx]])
        cmap = 'inferno'

        im1 = ax1.imshow(covs[mask_idx], vmin=min_val, vmax=max_val, cmap=cmap)
        ax1.set_xlabel(r'$\ell_1$')
        ax1.set_ylabel(r'$\ell_2$')
        ax1.set_title('Numerical covariance matrix')

        im2 = ax2.imshow(covs_th[mask_idx], vmin=min_val, vmax=max_val, cmap=cmap)
        ax2.set_xlabel(r'$\ell_1$')
        ax2.set_title('Theory covariance matrix')

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.05, 0.7])
        fig.colorbar(im2, label=r'$\textrm{Cov}[\ell_1, \ell_2]$', cax=cbar_ax)

        fig.suptitle(r'Comparison for mask {msk} ($f_\textrm{{sky}}: {fsky:.2f} \%$)'
                     .format(msk=mask_idx + 1, fsky=self.masks_f_sky[0] * 100))

        # * Plot the absolute difference

        print('Max difference {val:.3e}'.format(val=np.max(covs[mask_idx] - covs_th[mask_idx])))
        print('Min difference {val:.3e}'.format(val=np.min(covs[mask_idx] - covs_th[mask_idx])))

        fig3, ax3 = plt.subplots()
        im3 = ax3.imshow(covs[mask_idx] - covs_th[mask_idx], cmap=make_planck_colour_map(),
                         vmin=-np.max(covs[mask_idx] - covs_th[mask_idx]))
        ax3.set_xlabel(r'$\ell_1$')
        ax3.set_ylabel(r'$\ell_2$')
        ax3.set_title(f'Difference between covariance matrices for mask {mask_idx + 1}')
        fig3.colorbar(im3, label='Numerical - Theory')

        # * Plot the relative difference (only around diagonal)

        tmp_cov = covs[mask_idx]
        tmp_cov_th = covs_th[mask_idx]

        fig3, ax3 = plt.subplots()
        im3 = ax3.imshow(tmp_cov / tmp_cov_th - 1, cmap=make_planck_colour_map(),
                         vmin=-(np.max(tmp_cov / tmp_cov_th) - 1))
        ax3.set_xlabel(r'$\ell_1$')
        ax3.set_ylabel(r'$\ell_2$')
        ax3.set_title(f'Relative difference between covariance matrices for mask {mask_idx + 1}')
        fig3.colorbar(im3, label='Numerical / Theory - 1')

        # Now compute what the relative difference matrix is
        matrix = tmp_cov / tmp_cov_th - 1

        # Extract the diagonal and off-diagonals of this matrix
        diag = np.diagonal(matrix, offset=0)
        diag_m1 = np.diagonal(matrix, offset=-1)
        diag_m2 = np.diagonal(matrix, offset=-2)
        diag_p1 = np.diagonal(matrix, offset=+1)
        diag_p2 = np.diagonal(matrix, offset=+2)

        # Redefine ells here to start at zero
        ells = np.arange(0, lmax - 1)

        # New figure for the one diagonal and four off-diagonal components
        fig3, (ax1, ax2, ax3, ax4, ax5) = plt.subplots(figsize=(11, 6), nrows=5, sharex=True)
        ax1.scatter(ells[1:-1], diag_p2, c='lime', s=4)
        ax2.scatter(ells[1:], diag_p1, c='purple', s=4)
        ax3.scatter(ells, diag, c='cornflowerblue', s=4)
        ax4.scatter(ells[1:], diag_m1, c='orange', s=4)
        ax5.scatter(ells[2:-1], diag_m2[1:], c='hotpink', s=4)

        ax1.set_ylabel('Diag + 2')
        ax2.set_ylabel('Diag + 1')
        ax3.set_ylabel('Diag')
        ax4.set_ylabel('Diag - 1')
        ax5.set_ylabel('Diag - 2')

        ax5.set_xlabel(r'$\ell$')
        ax1.set_title('Numerical covariance / theory covariance - 1')
        fig3.tight_layout()

        # * Now plot how close to symmetric the covariance matrices are
        fig, ax = plt.subplots()

        min_val = np.min(covs_th[mask_idx] - covs_th[mask_idx].T)
        max_val = np.max(covs_th[mask_idx] - covs_th[mask_idx].T)
        cmap = make_planck_colour_map()

        im = ax.imshow(covs_th[mask_idx] - covs_th[mask_idx].T, vmin=min_val, vmax=max_val, cmap=cmap)

        ax.set_xlabel(r'$\ell_1$')
        ax.set_ylabel(r'$\ell_2$')
        ax.set_title('Anti-symmetry of the theory covariance matrix')
        fig.colorbar(im, label=r'$\textrm{Cov}[\ell_1, \ell_2] - \textrm{Cov}^\textrm{T}[\ell_1, \ell_2]$', aspect=15)
        fig.tight_layout()

        plt.show()

        fig, ax = plt.subplots(figsize=(11, 6))

        ax.loglog(ells, np.diagonal(covs[0]), c='tab:blue', lw=3, label='Numerical covariance')
        ax.loglog(ells, num_vars[0], c='orange', lw=1.5, label='Numerical variance', ls='--')
        ax.loglog(ells, np.diagonal(covs_th[0]), c='hotpink', lw=1.5, label='Theory covariance')

        # Also plot the expected variance of just the cosmic variance
        ax.loglog(ells, 2 / (2 * ells + 1) * (self.c_ells['W4xW4'][2:] * self.masks_f_sky[0]) ** 2, ls='--', lw=2,
                  color='cyan', label='Cosmic variance')

        ax.set_xlabel(r'$\ell$')
        ax.set_ylabel(r'$\textrm{Var}[C_\ell]$')
        ax.set_title(r'Comparison of $C_\ell$ variance as a function of $\ell$ between theory and numerics')

        plt.legend()
        plt.tight_layout()

        norm = mpl.colors.LogNorm(vmin=min(self.masks_f_sky), vmax=max(self.masks_f_sky))
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap='plasma')
        cmap.set_array([])

        fig, ax = plt.subplots(figsize=(11, 6))
        for idx in range(9):
            ax.semilogx(ells, np.diagonal(covs[idx]) / np.diagonal(covs_th[idx]) - 1,
                        c=cmap.to_rgba(self.masks_f_sky[idx]))

        ax.set_xlabel(r'$\ell$')
        ax.set_ylabel(r'$\textrm{Numerical Var}[C_\ell] / \textrm{Theory Var}[C_\ell] - 1$')
        ax.set_title('Ratio of the numetical to the theoretical variance')

        fig.colorbar(cmap, label=r'$f_\textrm{sky}$')
        plt.tight_layout()
        plt.show()

        # * Now make an animation of the theory covariance matrix
        fig, ax = plt.subplots()
        cmap = 'inferno'

        # Create the initial heatmap
        im = ax.imshow(covs[0], cmap=cmap)

        # Create a sub-axes that's on the right for the colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes('right', size='5%', pad=0.1)

        # fig.colorbar(im, label=r'$\textrm{Cov}[\ell_1, \ell_2]$', cax=cax)

        def animate(frame_num):
            """
                Function that gets called every time a new frame is drawn in the animation

            Args:
                frame_num (int): The current frame number, provided by matplotlib FuncAnimation function. Starts at 0

            Returns:
                None
            """
            # Clear the existing data on the axes
            ax.clear()
            cax.clear()

            # Update the heatmap
            im2 = ax.imshow(covs[frame_num], cmap=cmap)

            fig.colorbar(im2, label=r'$\textrm{Cov}[\ell_1, \ell_2]$', cax=cax)

            # Re-set axes labels
            ax.set_xlabel(r'$\ell_1$')
            ax.set_ylabel(r'$\ell_2$')

            # Also update the plot title with mask number and f_sky
            ax.set_title('Mask ' + str(frame_num + 1) + r', fsky: {fsky:.2f} \%'.format(
                fsky=self.masks_f_sky[frame_num] * 100))

        # Use the maptlotlib animation module to create the animation
        anim = pltanim.FuncAnimation(fig, animate, frames=len(covs), interval=750, repeat=True)

        plt.show()

    def plot_multiple_run_data(self, used_mask=False):
        """
        Function that plots data that has been calculated for multiple runs of Flask

        Args:
            used_mask (bool): Bool value to determine if a mask was used then computing the Cl coefficients,
            and so to plot the distribution of the masked Cl values too

        Returns:
            None
        """
        # Import the previously saved data
        data_df = pd.read_csv(self.folder_path + 'AggregateCls.csv')
        data_df_kde = pd.read_csv(self.folder_path + 'AggregateCls.csv')

        print('Total number of samples here is: ' + str(len(data_df)))

        # Lists where computed data will be stored into
        mean_cls = []
        var_cls = []

        # Skew and kurtosis computed explicitly from the raw data
        skew_cls = []
        kurt_cls = []

        # Skew and kurtosis where we numerically evaluate the moments, but divide by the expected variance
        skew_cls2 = []
        kurt_cls2 = []

        # Here, we want to find what the average and variance of the Cl's are at each ell
        for label, c_ells in data_df.items():
            # Calculate the mean and append it to our list
            mean_val = np.mean(c_ells)
            mean_cls.append(np.mean(c_ells))

            # Calculate the variance, normalise it through the mean value squared, and append it to our list
            var_cls.append(np.var(c_ells) / (mean_val * mean_val))

            # Also calculate the skew and kurtosis and append them to the lists too
            skew_cls.append(scistats.skew(np.array(c_ells), bias=True))
            kurt_cls.append(scistats.kurtosis(np.array(c_ells), bias=True))

            ell = int(label)
            exp_var = 2 / (2 * ell + 1)
            exp_var *= self.raw_c_ells['W2xW2'][ell] ** 2

            skew_cls2.append(scistats.moment(np.array(c_ells), moment=3) / (exp_var ** (3 / 2)))
            kurt_cls2.append((scistats.moment(np.array(c_ells), moment=4) / (exp_var ** 2)) - 3)

        # If want to plot masked data as well, then read it in now
        if used_mask:
            masked_data_df = pd.read_csv(self.folder_path + 'AggregateMaskedCls.csv')
            masked_data_df_kde = pd.read_csv(self.folder_path + 'AggregateMaskedCls.csv')

            print('Total number of samples in masked DataFrame is: ' + str(len(masked_data_df)))

            mean_cls_2 = []
            var_cls_2 = []

            skew_cls_2 = []
            kurt_cls_2 = []

            for label, c_ells in masked_data_df.items():
                # Calculate the mean and append it to our list
                mean_val = np.mean(c_ells)
                mean_cls_2.append(np.mean(c_ells))

                # Calculate the variance, normalise it through the mean value squared, and append it to our list
                var_cls_2.append(np.var(c_ells) / (mean_val * mean_val))

                # Also calculate the skew and kurtosis and append them to the lists too
                skew_cls_2.append(scistats.skew(np.array(c_ells), bias=True))
                kurt_cls_2.append(scistats.kurtosis(np.array(c_ells), bias=True))

        # Get the keys of the data frame
        keys = data_df.keys()

        # Dictionary which the Pearson data will be stored into
        pearson_data = {'x': [], 'y': [], 'r': []}

        covariance_data = {'x': [], 'y': [], 'raw_C': [], 'norm_C': []}

        # Go through each Cl combination and work out the Pearson correlation coefficient
        for ell1, c_ells1 in masked_data_df.items():
            # Only go up to ell1 of 125, otherwise figure too crowded
            if int(ell1) > 126:
                continue

            # Only want to compute samples where ell2 < ell1
            for ell2, c_ells2 in masked_data_df.items():
                if int(ell2) >= int(ell1):
                    continue

                # Use the Pearson test to get correlation coefficient
                r_val, p_val = scistats.pearsonr(c_ells1, c_ells2)

                # Store the calculation results in the Pearson data dictionary
                pearson_data['x'].append(int(ell1))
                pearson_data['y'].append(int(ell2))
                pearson_data['r'].append(r_val)

                # * Now manually compute the covariance matrix

                val = 0

                # The column in the c_ells dictionary that has lensing information stored in it
                camb_col = 'W4xW4' if self.galaxy_dens else 'W2xW2'

                # Go through each ell1 and ell2 list and compute the sum at that point
                for val1, val2 in zip(c_ells1.to_numpy(), c_ells2.to_numpy()):
                    val += (val1 - self.c_ells[camb_col][int(ell1)]) * \
                           (val2 - self.c_ells[camb_col][int(ell2)])

                # Normalise the sum to the number of data points
                val /= len(c_ells1)

                # Store the data in the covariance data dictionary
                covariance_data['x'].append(int(ell1))
                covariance_data['y'].append(int(ell2))
                covariance_data['raw_C'].append(val)

                # Now normalise through the intrinsic standard deviation at each ell value
                val /= self.c_ells[camb_col][int(ell1)] * np.sqrt(2 / (2 * int(ell1) + 1))
                val /= self.c_ells[camb_col][int(ell2)] * np.sqrt(2 / (2 * int(ell2) + 1))

                # Store that in the dictionary too
                covariance_data['norm_C'].append(val)

        # Convert the covariance dictionary to a dataframe
        covariance_df = pd.DataFrame(covariance_data)

        # Pivot the dataframe to get in correct form for plotting a heatmap
        raw_covariance_df = covariance_df.pivot('x', 'y', 'raw_C')
        norm_covariance_df = covariance_df.pivot('x', 'y', 'norm_C')

        # Turn our dictionary into a data-frame
        pearson_df = pd.DataFrame(pearson_data)

        # Pivot the data-frame to get it in the right format for plotting
        pearson_df = pearson_df.pivot('x', 'y', 'r')

        # Use a heatmap to plot the correlation coefficients
        with sns.axes_style("white"):
            plt.figure(figsize=(9, 8))
            ax = sns.heatmap(data=pearson_df, square=True, cmap="jet",
                             cbar_kws={'label': 'Pearson $r$ correlation coefficient'})
            ax.set_xlabel(r'$\ell_1$')
            ax.set_ylabel(r'$\ell_2$')
            plt.tight_layout()

            plt.figure(figsize=(9, 8))
            ax = sns.heatmap(data=raw_covariance_df, square=True, cmap="jet",
                             cbar_kws={'label': 'Covariance matrix $C_{ij}$'})
            ax.set_xlabel(r'$\ell_i$')
            ax.set_ylabel(r'$\ell_j$')
            ax.set_title('Raw covariance matrix')
            plt.tight_layout()

            plt.figure(figsize=(9, 8))
            ax = sns.heatmap(data=norm_covariance_df, square=True, cmap="seismic", center=0,
                             cbar_kws={'label': r'Covariance matrix $C_{ij} / \sigma_i \sigma_j$'})
            ax.set_xlabel(r'$\ell_i$')
            ax.set_ylabel(r'$\ell_j$')
            ax.set_title('Normalised covariance matrix for masked convergence')
            plt.tight_layout()

            plt.show()

        # In order to plot a KDE, we need to increase the size of the values in order to make them order-unity
        for idx in range(len(keys)):
            data_df_kde[data_df_kde.keys()[idx]] = (data_df_kde[data_df_kde.keys()[idx]] * 1E7)

        grid2 = sns.PairGrid(data_df_kde, vars=[keys[0], keys[2], keys[5], keys[10], keys[15], keys[25], keys[50]])
        grid2.map_diag(sns.histplot)
        grid2.map_lower(sns.kdeplot, shade=True, levels=4)
        grid2.map_upper(sns.histplot)
        grid2.tight_layout()
        plt.show(block=False)

        # Also plot the masked KDE plot, if we used one
        if used_mask:
            # Also make all values order unity by multiplying
            for idx in range(len(keys)):
                masked_data_df_kde[masked_data_df_kde.keys()[idx]] = \
                    (masked_data_df_kde[masked_data_df_kde.keys()[idx]] * 1E7)

            grid2 = sns.PairGrid(masked_data_df_kde,
                                 vars=[keys[0], keys[2], keys[5], keys[10], keys[15], keys[25], keys[50]])
            grid2.map_diag(sns.histplot)
            grid2.map_lower(sns.kdeplot, shade=True, levels=4)
            grid2.map_upper(sns.histplot)
            grid2.tight_layout()
            plt.show(block=False)

        # * Print the raw data for each Cl dataset - very messy but gets the point across
        plt.figure(figsize=(13, 7))
        for cl in data_df.itertuples(index=False):
            plt.loglog(self.ells, cl, alpha=0.5, linewidth=0.75)

        plt.loglog(self.ells, mean_cls, lw=2, color='tab:blue', label=r'Average $C_\ell$')
        if used_mask:
            plt.loglog(self.ells, mean_cls_2, lw=2, color='hotpink', label=r'Average masked $C_\ell$')
        plt.loglog(self.ells, self.c_ells['W4xW4'][2:], lw=1, color='tab:cyan', label=r'CAMB $C_\ell$')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')

        plt.legend()
        plt.tight_layout()

        # * Plot the deviation of the mean_cls vs input cls from CAMB
        camb_col = 'W4xW4' if self.galaxy_dens else 'W2xW2'

        plt.figure(figsize=(13, 7))
        plt.semilogx(self.ells, (np.array(mean_cls) / self.c_ells[camb_col][2:]) - 1, lw=2, color='tab:blue',
                     label=r'Unmasked $C_\ell$')
        if used_mask:
            plt.semilogx(self.ells, (np.array(mean_cls_2) / self.c_ells[camb_col][2:]) - 1, lw=2, color='hotpink',
                         label=r'Masked $C_\ell$')

        plt.title(r'Deviation of the average $C_\ell$ with respect to the input $C_\ell$ values')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$ (C_\ell^\textrm{Avg} - C_\ell^\textrm{In} ) / C_\ell^\textrm{In} $')

        plt.legend()
        plt.tight_layout()

        # Also plot the skew
        plt.figure(figsize=(13, 7))

        # Plot the skew
        plt.semilogx(self.ells, skew_cls, lw=1, color='tab:blue', label=r'Full $C_\ell$')
        if used_mask:
            plt.semilogx(self.ells, skew_cls_2, lw=1, color='hotpink', label=r'Masked $C_\ell$')

        # Plot the expected skew for a Gamma distribution
        k = (2 * self.ells + 1) / 2
        plt.semilogx(self.ells, 2 / np.sqrt(k), color='purple', lw=2, label=r'$\Gamma$ function prediction')

        plt.title('Skew')
        plt.xlabel(r'$\ell$')
        plt.legend()
        plt.tight_layout()

        # And kurtosis
        plt.figure(figsize=(13, 7))
        plt.semilogx(self.ells, kurt_cls, lw=1, color='tab:blue', label=r'Full $C_\ell$')
        if used_mask:
            plt.semilogx(self.ells, kurt_cls_2, lw=1, color='hotpink', label=r'Masked $C_\ell$')

        plt.semilogx(self.ells, 6 / k, color='purple', lw=2, label=r'$\Gamma$ function prediction')
        plt.title('Kurtosis')
        plt.xlabel(r'$\ell$')
        plt.legend()
        plt.tight_layout()

        # Now plot the variance of the Cl's with the expected values from the Cosmic variance
        plt.figure(figsize=(13, 7))
        plt.loglog(self.ells, var_cls, lw=1, color='tab:blue', label=r'Full $C_\ell$')
        if used_mask:
            plt.loglog(self.ells, var_cls_2, lw=1, color='hotpink', label=r'Masked $C_\ell$')

        plt.loglog(self.ells, 2 / (2 * self.ells + 1), lw=2, color='purple', label='Cosmic variance')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\textrm{Var}[C_\ell] / \textrm{Avg}[C_\ell]^2$')
        plt.title('Variance divided by the squared-average for the masked and unmasked convergence field')

        plt.legend()
        plt.tight_layout()
        plt.show()

        return

        # * Code that computes the ks-test to see if our Cl distributions *are* statistically Gaussian or Gamma-func

        for ell, c_ells in data_df.items():
            stdev = self.raw_c_ells['W2xW2'][int(ell)] * np.sqrt(2 / (2 * int(ell) + 1))
            k_stat, p_val = scistats.kstest(c_ells, 'norm', args=[self.raw_c_ells['W2xW2'][int(ell)], stdev])
            print(k_stat, p_val)

            k_stat, p_val = scistats.kstest(c_ells, 'gamma', args=[(2 * int(ell) + 1) / 2,
                                                                   self.raw_c_ells['W2xW2'][int(ell)] *
                                                                   2 / (2 * int(ell) + 1)])
            print(k_stat, p_val)

            plt.figure(figsize=(13, 7))
            plt.hist(c_ells, bins=50)

            norm_vals = np.random.normal(loc=self.raw_c_ells['W2xW2'][int(ell)], scale=stdev, size=len(c_ells))
            plt.hist(norm_vals[norm_vals > 0], bins=50, alpha=0.5)

            gamma_vals = np.random.gamma(shape=(2 * int(ell) + 1) / 2,
                                         scale=self.raw_c_ells['W2xW2'][int(ell)] * 2 / (2 * int(ell) + 1),
                                         size=len(c_ells))
            plt.hist(gamma_vals, bins=50, alpha=0.5)

            plt.title('l = ' + str(ell))
            plt.show()

    def plot_ridge_plot(self):
        """
        Function to plot a ridge-plot of the Cl values

        Returns:
            None
        """

        # Read in the data to a Pandas DataFrame
        data_df = pd.read_csv(self.folder_path + 'AggregateCls.csv')
        masked_data_df = pd.read_csv(self.folder_path + 'AggregateMaskedCls.csv')

        # Normalise each value by 1E7 to get working KDE plots
        for idx in range(len(data_df.keys())):
            data_df[data_df.keys()[idx]] = (data_df[data_df.keys()[idx]] * 1E7)

        for idx in range(len(masked_data_df.keys())):
            masked_data_df[masked_data_df.keys()[idx]] = (masked_data_df[masked_data_df.keys()[idx]] * 1E7)

        data_frames = []

        # Go through each ell and extract the Cl data
        for label, c_ells in data_df.items():
            # Subject to the condition where ell <= 25
            if int(label) > 25:
                continue
            tmp_df = pd.DataFrame(
                    {'ell': label, 'Cl': (c_ells - np.mean(c_ells)) / np.mean(c_ells), 'Type': 'Unmasked'})
            data_frames.append(tmp_df)

        for label, c_ells in masked_data_df.items():
            # Subject to the condition where ell <= 25
            if int(label) > 25:
                continue
            tmp_df = pd.DataFrame({'ell': label, 'Cl': (c_ells - np.mean(c_ells)) / np.mean(c_ells), 'Type': 'Masked'})
            data_frames.append(tmp_df)

        data = pd.concat(data_frames, ignore_index=True)

        # Resets the Seaborn pallet to have a white figure background
        sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

        # Create the colour map for our ell values
        # pal = sns.cubehelix_palette(len(data_frames), rot=-0.2, light=0.75, dark=0.05, gamma=1.25)
        pal = sns.cubehelix_palette(2, rot=-0.2, light=0.85, dark=0.25, gamma=1.25)

        # Initialize the FacetGrid object
        g = sns.FacetGrid(data, row="ell", hue="Type", aspect=20, height=0.5, palette=pal)

        # Draw the densities in a few steps
        g.map(sns.kdeplot, "Cl", bw_adjust=0.5, clip_on=False, fill=True, alpha=0.6, linewidth=2)
        # ? g.map(sns.kdeplot, "Cl", clip_on=False, color="w", lw=2, bw_adjust=0.5)
        g.map(plt.axhline, y=0, lw=2, clip_on=False)
        # g.map(plt.axhline, y=0, xmin=-0.75, xmax=0.75,  lw=2, clip_on=False)
        g.map(plt.axvline, x=0, linestyle='--', color='w', lw=1, alpha=0.8, clip_on=False)

        # Define and use a simple function to label the plot in axes coordinates
        def label_func(x, color, label):
            ax = plt.gca()

            global num_ell
            num_ell += 0.5
            label = num_ell

            if label == int(label):
                ax.text(0, .2, int(label), fontweight="bold", color=color,
                        ha="left", va="center", transform=ax.transAxes)

        g.map(label_func, "Cl")

        # Function to plot the expected values for the Γ-distribution which describes the distribution of the Cl's
        def plot_gamma(x, color, label):
            # Get the current plot's axis
            ax = plt.gca()

            # Compute the expected variance from the ell value
            var = 2 / (2 * int(label) + 1)

            # Sort the x values ascending so our plot is smooth and continuous
            x = np.sort(x)

            # Plot the Γ-distribution PDF given the variance for the ell value, in blue.
            ax.plot(x, scistats.gamma.pdf(x, a=1 / var, loc=-1, scale=var), lw=3, color='tab:blue')

            ax.plot(x, scistats.norm.pdf(x, loc=0, scale=np.sqrt(var)), lw=2, color='tab:green')

        # Map the above function to the Cl data
        # g.map(plot_gamma, "Cl")

        # Set the subplots to overlap
        g.fig.subplots_adjust(hspace=-0.45)

        # Add title to figure and adjust its position
        g.fig.subplots_adjust(top=0.995)
        g.fig.suptitle(r'Distribution of the $C_\ell$ values', fontsize=18)
        plt.xlabel(r'$[C_\ell - \bar{C}_\ell] / \bar{C}_\ell$', fontsize=16)
        g.set(xlim=(-1, 1))

        # Remove axes details that don't play well with overlap
        g.set_titles("")
        g.set(yticks=[])
        g.despine(bottom=True, left=True)

        g.add_legend(fontsize='large')
        plt.show()

    @staticmethod
    def get_file_num_lines(input_file):
        """
        Function to get the number of lines for a given input file

        Args:
            input_file (str): Path to input file

        Returns:
            (int) The number of lines
        """

        # Ensure that the file exists before opening it
        if not os.path.isfile(input_file):
            raise RuntimeError('Cannot find the given file path! Please check that the input is correct')

        # Open the file in a read-only format
        with open(input_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                pass

        # Return
        return line_num

    def trim_flask_alm_output(self):
        """
        By default, Flask outputs the Alm file with a header.
        However, for this file to be read in by the C++ code we need to trim off this header.

        :return: None
        """

        # Go through each redshift bin
        for z in range(1, self.num_redshift_bins + 1):
            # Get the path for the input file, and a dummy temporary file
            input_file = self.folder_path + 'Output-shear-alm-f1z' + str(z) + '.dat'
            dummy_file = self.folder_path + 'dummy_file.dat'

            # Open the input file to get the number of lines in the file
            line_num = self.get_file_num_lines(input_file)

            reject_lines = [0, 1]

            if line_num == hp.sphtfunc.Alm.getsize(self.ell_max) - 1:
                with open(input_file, 'r') as read_obj, open(dummy_file, 'w') as write_obj:
                    # Line by line copy data from original file to dummy file
                    for line_num, line in enumerate(read_obj):
                        if line_num not in reject_lines:
                            write_obj.write(line)

                os.remove(input_file)
                os.rename(dummy_file, input_file)

    def estimate_cl_from_alm(self):
        """
        Function that reads in the alm values as outputted by Flask and converts them to Cl's

        :return: None, as we store the Cl information in the "self.my_cls" class property
        """

        # Go through each z bin
        for z in range(1, self.num_redshift_bins + 1):
            start_time = time.time()

            input_file = self.folder_path + 'Output-shear-alm-f1z' + str(z) + '.dat'

            # Read in the alm values for the z bin into a Pandas DataFrame
            alms = pd.read_csv(input_file, sep=r'\s+',
                               names=['l', 'm', 'Re', 'Im'],
                               header=None,
                               skiprows=2 if
                               self.get_file_num_lines(input_file) == hp.sphtfunc.Alm.getsize(self.ell_max) - 1
                               else 0,
                               dtype={'l': np.uint16, 'm': np.uint16, 'Re': np.float64, 'Im': np.float64})

            # Initialise an array for storing the Cl's in, default to zero
            cls = np.zeros(self.ell_max + 1)

            # TODO: check that this Cl calculation is consistent with Flask's notation
            # I think they only divide by ell + 1, and don't bother with a factor of two.

            # Go through each row in the alm dataframe
            for row in alms.itertuples(index=False):
                # Compute the modulus-square of the alm coefficient
                alm_sq = row.Re * row.Re + row.Im * row.Im

                # If this alm has m != 0, then we can multiply by two from symmetries of the alm
                cls[row.l] += (alm_sq if row.m == 0 else 2 * alm_sq)

            # Normalise our average by dividing by (2 * ell + 1)
            for ell in self.ells:
                cls[ell] /= 2 * ell + 1

            finish_time = time.time()
            print('My Cl calculation took {num:.2f} seconds'.format(num=finish_time - start_time))

            # Store our Cl list in the class
            self.my_cls.append(cls)

    def plot_map_to_alm_diff(self):
        maps = hp.read_map(self.folder_path + 'Output-shear-map-fits-f1z1.fits', verbose=False, field=None)

        # Split the shae map into a convergence, shear1, and shear2 maps
        converg_map = maps[0]
        shear_map1 = maps[1]
        shear_map2 = maps[2]

        shear_cl1 = hp.sphtfunc.anafast(shear_map1, lmax=self.ell_max)

        plt.figure(figsize=(12, 7))
        plt.loglog(self.ells, self.ells * (self.ells + 1) * self.my_cls[0][2:] / (2 * np.pi),
                   label='My Cls', color='navy', lw=1.75)
        plt.loglog(self.ells, self.ells * (self.ells + 1) * shear_cl1[2:] / (2 * np.pi),
                   label=r'Cls from Flask', color='steelblue', lw=0.5, alpha=0.5)
        plt.title(r'Comparison of $C_\ell$ computation for redshift bin 1')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        plt.legend()
        plt.xlim(right=1E3)
        plt.tight_layout()
        plt.show()

    def plot_flask_output(self):
        """
        This function provides utilities to read in data produced by Flask and plot then using HealPy to plot the
        output maps and power spectra from maps

        Returns:
            None
        """

        # Obtain the Planck colour-map used for plotting maps here
        planck_cmap = make_planck_colour_map()

        # If we have galaxy counts, then the lensing signal is the second field, else it is simply the first
        field_num = 2 if self.galaxy_dens else 1

        # Find the header character of the output Cl file from Flask
        cl_head_char = subprocess.run('head -c 1 Output-Cl.dat',
                                      stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                                      cwd=self.folder_path, shell=True).stdout

        if cl_head_char == '#':
            subprocess.run('tail -c +3 Output-Cl.dat >> tmp_file.dat && mv tmp_file.dat Output-Cl.dat',
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True,
                           cwd=self.folder_path, shell=True)

        flask_cl_df = pd.read_csv(self.folder_path + 'Output-Cl.dat', sep=r'\s+')
        flask_cl_df.rename(columns={"l": "ell"}, inplace=True)

        # Go through each redshift bin individually and make plots for each
        for z_bin in range(1, self.num_redshift_bins + 1):
            # Read in the generated shear map with noise from Flask
            converg_w_noise = hp.read_map(
                    self.folder_path + 'Output-poisson-map-f' + str(field_num) + 'z' + str(z_bin) + '.fits',
                    verbose=False, field=None)

            ellip_w_noise = hp.read_map(
                    self.folder_path + 'Output-Ellip-map-f' + str(1) + 'z' + str(z_bin) + '.fits',
                    verbose=False, field=None)

            # Split the shae map into a convergence, shear1, and shear2 maps
            # converg_map = maps[0]
            # shear_map1 = maps[1]
            # shear_map2 = maps[2]

            # Plot the convergence map
            hp.mollview(converg_w_noise, cmap=planck_cmap,
                        title='Convergence with shape noise for redshift bin ' + str(z_bin))
            hp.graticule(verbose=False, alpha=0.6)

            hp.mollview(ellip_w_noise[1], cmap=planck_cmap, title='Ellipticity map for redshift bin ' + str(z_bin))
            hp.graticule(verbose=False, alpha=0.6)

            # Use HealPy functions to turn the maps into Cl's
            start_time = time.time()
            converg_cl = hp.sphtfunc.anafast(converg_w_noise, lmax=self.ell_max)
            print('Cl estimation fromt the convergence map w/noise took {sec} seconds'
                  .format(sec=time.time() - start_time))

            # Now compute the ellipticity Cl's
            ellip_cl1 = hp.sphtfunc.anafast(ellip_w_noise[1], lmax=self.ell_max)
            ellip_cl2 = hp.sphtfunc.anafast(ellip_w_noise[2], lmax=self.ell_max)

            """
            start_time = time.time()
            # ? shear_cl1 = hp.sphtfunc.anafast(shear_map1, lmax=self.ell_max)
            print('Shear1 took {sec} seconds'.format(sec=time.time() - start_time))

            start_time = time.time()
            # ? shear_cl2 = hp.sphtfunc.anafast(shear_map2, lmax=self.ell_max)
            print('Shear2 took {sec} seconds'.format(sec=time.time() - start_time))
            """

            # Plot various Cl's
            plt.figure(figsize=(13, 7))

            # Plot the theory Cl's
            plt.loglog(self.ells, self.ells * (self.ells + 1) *
                       self.raw_c_ells['W' + str(field_num * z_bin) + 'x' + 'W' + str(field_num * z_bin)][2:]
                       / (2 * np.pi),
                       label=r'$C_\ell$ input', lw=1.5, color='navy')

            # Plot the full Cl's recovered from the map
            plt.loglog(self.ells, self.ells * (self.ells + 1) * converg_cl[2:] / (2 * np.pi),
                       label=r'Converg w/shp. noise', color='tab:green')

            # Plot the Cl's that just have cosmic variance, no shape noise
            plt.loglog(flask_cl_df['ell'],
                       flask_cl_df['ell'] * (flask_cl_df['ell'] + 1) *
                       flask_cl_df['Cl-f' + str(field_num) + 'z' + str(z_bin) + 'f' + str(field_num) + 'z' + str(z_bin)]
                       / (2 * np.pi), lw=2,
                       color='tab:blue', label='Converg w/out shp. noise')

            plt.loglog(self.ells, self.ells * (self.ells + 1) * ellip_cl1[2:] / (2 * np.pi),
                       label=r'Ellip1', color='yellow')

            plt.loglog(self.ells, self.ells * (self.ells + 1) * ellip_cl2[2:] / (2 * np.pi),
                       label=r'Ellip2', color='hotpink')

            plt.title(r'$C_\ell$ for redshift bin ' + str(z_bin))
            plt.xlabel(r'$\ell$')
            plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
            plt.legend()
            plt.tight_layout()
            # plt.show()

            # Plot various differences in the Cl's
            plt.figure(figsize=(13, 7))

            # Plot the difference between Cl's with shape noise and cosmic variance and input
            plt.loglog(self.ells, np.abs(converg_cl[2:] -
                       self.raw_c_ells['W' + str(field_num * z_bin) + 'x' + 'W' + str(field_num * z_bin)][2:]) /
                       self.raw_c_ells['W' + str(field_num * z_bin) + 'x' + 'W' + str(field_num * z_bin)][2:],
                       label='$a=$ w/shp. noise, $b=$ input', color='tab:green', lw=1.5)

            # Plot the difference between Cl's without shape noise (but cosmic variance) and input
            plt.loglog(self.ells, np.abs(
                    flask_cl_df['Cl-f' + str(field_num) + 'z' + str(z_bin) + 'f' + str(field_num) + 'z' + str(z_bin)] -
                    self.raw_c_ells['W' + str(field_num * z_bin) + 'x' + 'W' + str(field_num * z_bin)][2:]) /
                    self.raw_c_ells['W' + str(field_num * z_bin) + 'x' + 'W' + str(field_num * z_bin)][2:],
                    label='$a=$ w/out shp. noise, $b=$ input', color='tab:blue', lw=1.5, alpha=0.55)

            plt.xlabel(r'$\ell$')
            plt.ylabel(r'$|C_{\ell}^{a} - C_{\ell}^{b} | / C_{\ell}^{b} $')
            plt.title(r'Relative difference in $C_\ell$ for redshift bin ' + str(z_bin))
            plt.legend()
            plt.tight_layout()

        plt.show()

    def use_cpp_thingy(self):
        # Load the C++ shared library file that has already been compiled.
        lib = ctypes.CDLL('/home/amaraio/Documents/PhD/Codes/LensingMapMaking/lib/cpp/thingy.so')

        # Define a custom type: a pointer to a double
        c_double_p = ctypes.POINTER(ctypes.c_double)

        # Sets the arguments of the alm_to_clm function as a char pointer and double pointer
        lib.alm_to_cl.argtypes = [ctypes.c_char_p, c_double_p]

        # Get the size of the cl array
        size = self.ell_max + 1

        for z in range(1, self.num_redshift_bins + 1):
            # Initialise a Cl array of double pointers with length size
            cls = (ctypes.c_double * size)()

            # Set up the path variable that gets passed to the C++ code for reading in the data
            # Note that we have to use a "bytes" string for it to be read correctly
            file_path = ctypes.c_char_p(
                    b'/home/amaraio/Documents/PhD/Codes/LensingMapMaking/Data/Non-linear/Output-shear-alm-f1z'
                    + str(z).encode() + b'.dat')

            start_time = time.time()

            # Hand off to the C++ code
            lib.alm_to_cl(file_path, ctypes.cast(cls, c_double_p))

            # See how long it took
            print('My cpp Cl calculation took {num:.2f} seconds'.format(num=time.time() - start_time))

            # Store the cpp Cl's in the class
            self.my_cls_cpp.append(cls)

    @staticmethod
    def use_cpp_map_to_alm():
        # Load the C++ shared library file that has already been compiled.
        lib = ctypes.CDLL('/home/amaraio/Documents/PhD/Codes/LensingMapMaking/lib/cpp/thingy.so')

        # Sets the arguments of the map_to_alm function as a char pointer
        lib.alm_to_cl.argtypes = [ctypes.c_char_p]

        # Set up the path variable that gets passed to the C++ code for reading in the data
        # Note that we have to use a "bytes" string for it to be read correctly
        file_path = ctypes.c_char_p(
                b'/home/amaraio/Documents/PhD/Codes/LensingMapMaking/Data/Non-linear/Output-shear-map.dat')

        start_time = time.time()

        # Hand off to the C++ code
        lib.map_to_alm(file_path)

        # See how long it took
        print('My cpp map to alm calculation took {num:.2f} seconds'.format(num=time.time() - start_time))

    @staticmethod
    def experimenting_with_masks():
        """
        Here, we are playing around with masks and seeing what effects they have on the data!
        Not a real science function

        Returns:
            None
        """

        wmap_data = hp.read_map('./resources/existing_maps/wmap_band_iqumap_r9_7yr_Q_v4.fits')
        wmap_mask = hp.read_map('./resources/existing_maps/wmap_temperature_analysis_mask_r9_7yr_v4.fits').astype(np.bool)

        wmap_masked = hp.ma(wmap_data) * 1E3
        wmap_masked.mask = np.logical_not(wmap_mask)

        wmap_masked += np.random.normal(loc=0, scale=250, size=12 * 512**2)

        planck_cmap = make_planck_colour_map()
        hp.mollview(wmap_data, cmap='jet', norm="hist", min=-1, max=1, title='WMAP data with no mask (hist. equal cmap)')
        hp.graticule(verbose=False, alpha=0.8)

        # Plot the WMAP mask using the two coordinate systems
        hp.mollview(wmap_mask, cmap='magma', title='WMAP mask in Galactic coordinates', coord='G')
        hp.graticule(verbose=False, alpha=0.8)

        hp.mollview(wmap_mask, cmap='magma', title='WMAP mask in Ecliptic coordinates', coord='GE')
        hp.graticule(verbose=False, alpha=0.8)

        hp.mollview(wmap_masked.filled(), cmap='viridis', title='WMAP data with mask', unit=r'$\mu$K')
        hp.graticule(verbose=False, alpha=0.8)
        plt.show()

    def euclid_masks(self):
        """
        Function that uses the preliminary Euclid masks and recovers a power spectrum using it

        Returns:
            None
        """

        # Read in the preliminary Euclid mask using the correct N_side
        euclid_mask = hp.read_map('./resources/Euclid_masks/Euclid-gal-mask-2048.fits', verbose=False).astype(
            np.bool)

        # Compute the fraction of sky let through by the mask
        sky_fraction = euclid_mask.sum() / euclid_mask.size

        print('Fraction of sky allowed through by the mask is {num:.2f} %'.format(num=100 * sky_fraction))

        # Plot the Euclid mask in galactic coordinates
        hp.mollview(euclid_mask, cmap='seismic', title='Euclid mask in Galactic coordinates', coord='G')
        hp.graticule(verbose=False, alpha=0.8)

        # Plot the Euclid mask in ecliptic coordinates
        hp.mollview(euclid_mask, cmap='seismic', title='Euclid mask in Ecliptic coordinates', coord='GE')
        hp.graticule(verbose=False, alpha=0.8)

        # Also plot the mask using a spherical model
        hp.orthview(euclid_mask, rot=(0, 20, 15), title="Euclid mask projected on a sphere", cmap='seismic')
        hp.graticule(verbose=False)

        # Read in the convergence map that has been previously calculated
        converg_map = hp.read_map(self.folder_path + 'Output-poisson-map-f2z2.fits', verbose=False, field=None)

        # Use the mask to convert the raw map into a masked map
        masked_map = hp.ma(converg_map)
        masked_map.mask = np.logical_not(euclid_mask)

        # Plot the raw and masked maps
        hp.mollview(converg_map, cmap=make_planck_colour_map(), title='Raw converg data')
        hp.graticule(verbose=False, alpha=0.8)

        hp.mollview(masked_map, cmap=make_planck_colour_map(), title='Converg data with mask')
        hp.graticule(verbose=False, alpha=0.8)

        # Convert the raw map to a set of Cl values
        converg_cls = hp.anafast(converg_map, lmax=self.ell_max)

        # Use our masked map to generate a set of Cl values
        masked_cls = hp.anafast(masked_map, lmax=self.ell_max)

        # * We now want to perform a spherical harmonic expansion of the mask, to find what Cl values it brings
        mask_cls = hp.anafast(euclid_mask, lmax=self.ell_max)

        ell = np.arange(2, self.ell_max + 1)

        cl_df = pd.read_csv(self.folder_path + 'Unmasked_Cls.dat', sep=r'\s+')
        cl_df_mask = pd.read_csv(self.folder_path + 'Masked_Cls.dat', sep=r'\s+')

        # Plot the recovered Cl values
        plt.figure(figsize=(13, 7))
        plt.loglog(ell, ell * (ell + 1) * converg_cls[2:] / (2 * np.pi),
                   lw=3, c='tab:blue', label=r'$C_\ell$ recovered from full map')
        plt.loglog(ell, ell * (ell + 1) * masked_cls[2:] / (2 * np.pi),
                   lw=3, c='tab:green', label=r'$C_\ell$ recovered from masked map')
        plt.loglog(ell, ell * (ell + 1) * masked_cls[2:] / (2 * np.pi) / sky_fraction,
                   lw=1.5, c='orange', label=r'Naïve correction')
        plt.loglog(ell, ell * (ell + 1) * self.raw_c_ells['W4xW4'][2:] / (2 * np.pi),
                   lw=2, c='tab:cyan', label=r'Input $C_\ell$')
        plt.loglog(cl_df['ell'], cl_df['ell'] * (cl_df['ell'] + 1) * cl_df['Cl-f2z2f2z2'] / (2 * np.pi),
                   lw=2, c='hotpink', label=r'Flask\'s $C_\ell$')
        plt.loglog(cl_df_mask['ell'],
                   cl_df_mask['ell'] * (cl_df_mask['ell'] + 1) * cl_df_mask['Cl-f2z2f2z2'] / (2 * np.pi),
                   lw=2, c='navy', label=r'Flask\'s $C_\ell$')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_{\ell}^{\kappa\kappa} / 2 \pi$')
        plt.title('Convergence power spectrum with and without mask')

        plt.legend()
        plt.tight_layout()

        plt.figure(figsize=(13, 7))
        plt.semilogx(ell, (masked_cls[2:] / sky_fraction) / converg_cls[2:],
                     lw=3, c='tab:blue', label=r'$C_\ell$ recovered from full map')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\hat{C}_{\ell}^{\kappa\kappa} / C_{\ell}^{\kappa\kappa} $')
        plt.title(r'Ratio of recovered $C_{\ell}^{\kappa\kappa}$ from mask to unmasked values')

        plt.tight_layout()

        # Now plot the Cls from the Euclid mask
        plt.figure(figsize=(13, 7))
        plt.loglog(ell, ell * (ell + 1) * mask_cls[2:] / (2 * np.pi),
                   lw=3, c='tab:blue', label=r'$C_\ell$ recovered from full map')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_{\ell} / 2 \pi$')
        plt.title(r'Power spectrum of the mask')

        plt.tight_layout()
        plt.show()

    def generate_set_of_masks(self, apodizing_scale=None):
        """
        Here, we want to generate a set of masks that have different f_sky values, and so their effects on the recovered
        power spectra will be different.

        Args:
            apodizing_scale (float): If set, determines the scale that the generated mask should be apodized over (in
                                     degrees), if None then do not apodize.

        Returns:
            None
        """
        # Turn the class's N_side parameter into the number of pixels in the map
        n_pix = 12 * self.n_side ** 2

        if not os.path.isdir(self.folder_path + 'Masks'):
            os.makedirs(self.folder_path + 'Masks')

        if not os.path.isdir(self.folder_path + 'MaskedOutput'):
            os.makedirs(self.folder_path + 'MaskedOutput')

        plt.figure(figsize=(17, 10.25))

        # Go through a range of theta values that constructs us the galactic and ecliptic plane cuts
        for idx, theta in enumerate([35, 55, 75, 95, 125, 175, 275, 405, 505], start=1):
            # Create the galactic plane cut coordinates
            gal_th1 = (theta / 2 - 0.5) * np.pi / theta
            gal_th2 = (theta / 2 + 0.5) * np.pi / theta

            theta += 2

            # Now create the ecliptic plane cut coordinates
            elp_th1 = (theta / 2 - 0.5) * np.pi / theta
            elp_th2 = (theta / 2 + 0.5) * np.pi / theta

            # Create cut in galactic plane
            region_gal = hp.query_strip(nside=self.n_side, theta1=gal_th1, theta2=gal_th2)
            # Create cut in ecliptic plane
            region_elp = hp.query_strip(nside=self.n_side, theta1=elp_th1, theta2=elp_th2)

            map_gal = np.ones(n_pix, dtype=np.bool)
            map_elp = np.ones(n_pix, dtype=np.bool)

            # Mask out any regions according to our mask
            map_gal[region_gal] = 0
            map_elp[region_elp] = 0

            # Combine both masks into a single map in galactic coordinates
            map_both = np.logical_and(map_gal, hp.rotator.Rotator(coord='EG').rotate_map_pixel(map_elp), dtype=np.bool)

            # Manually delete maps to manage memory
            del map_gal
            del map_elp

            # Compute the sky fraction allowed through by this mask, before apodizing
            f_sky = map_both.sum() / map_both.size

            if apodizing_scale is not None:
                print('Fraction of sky allowed through by the dummy mask {msk} before apodizing is: {num:.9f} %'.format(
                    msk=idx, num=f_sky * 100))

                # Import NaMaster to do the apodizing
                import pymaster as mst

                # Apodize the mask using the provided scale. Defaults to using type "C2"
                map_both = mst.utils.mask_apodization(map_both, aposize=apodizing_scale, apotype='C2')

                # Compute the sky fraction again, after apodizing
                f_sky = map_both.sum() / map_both.size

                print('Fraction of sky allowed through by the dummy mask {msk} after apodizing is: {num:.9f} %'.format(
                    msk=idx, num=f_sky * 100))

            else:
                print('Fraction of sky allowed through by the dummy mask {msk} is {num:.9f} %'.format(
                    msk=idx, num=f_sky * 100))

            # Store this in the class
            self.masks_f_sky.append(f_sky)

            # Plot the current mask onto the existing figure
            hp.mollview(map_both, title=r"Mask with fsky {num:.2f} \%".format(num=f_sky), cmap='seismic',
                        sub=[3, 3, idx], margins=(0.005, -0.5, 0.0, -0.5), xsize=400, cbar=False)
            hp.graticule(verbose=False)

            # Save the mask to the masks sub-folder
            hp.write_map(self.folder_path + 'Masks/Mask' + str(idx) + ('_apo' if apodizing_scale is not None else '') +
                         '.fits', map_both, overwrite=True)

        plt.show()

        # Append the Euclid mask's f_sky
        self.masks_f_sky.append(self.mask_f_sky)

        # Also want to use the unmasked sky too
        self.masks_f_sky.append(1)

    def generate_small_fsky_masks(self):
        """
        Function that generates a set of masks that are similar to the above function, but exclusively makes masks
        with only a small f_sky

        Returns:
            None
        """

        # Convert the N_side parameter to the number of pixels
        n_pix = 12 * self.n_side ** 2

        # Ensure that the Masks sub-directory exists, if not then make it
        if not os.path.isdir(self.folder_path + 'Masks'):
            os.makedirs(self.folder_path + 'Masks')

        if not os.path.isdir(self.folder_path + 'MaskedOutput'):
            os.makedirs(self.folder_path + 'MaskedOutput')

        # We want to only allow data through in a small central region around the centre of the mask
        # This converts the (theta, phi) values to a HealPy vector
        vec = hp.ang2vec(np.pi / 2, 0)

        plt.figure(figsize=(15, 7.25))

        for idx, angle in enumerate([0.5, 0.35, 0.2, 0.1, 0.075, 0.05, 0.025, 0.01, 0.005], start=1):
            # Initalise blank mask as an array of zeros
            mask = np.zeros(n_pix)

            # We now find the region around the centre of the map that subtended by the given angle
            region = hp.query_disc(self.n_side, vec, angle)

            # Allow data to pass through in this region
            mask[region] = 1

            # Compute the f_sky parameter
            f_sky = mask.sum() / mask.size

            # Store this in the class
            self.masks_f_sky.append(f_sky)

            print('Fraction of sky allowed through by the dummy mask {msk} is {num:.9f} %'.format(msk=idx,
                                                                                                  num=f_sky * 100))

            # Plot this mask on a grid of masks
            hp.mollview(mask, title=r"Mask with fsky {num:.4f} \%".format(num=f_sky * 100), cmap='seismic',
                        sub=[3, 3, idx], margins=(0.005, -0.5, 0.0, -0.5), xsize=400, cbar=False)
            hp.graticule(verbose=False)

            # Save this mask to the folder
            hp.write_map(self.folder_path + 'Masks/Mask' + str(idx) + '.fits', mask, overwrite=True)

        plt.show()

        # Append the Euclid mask's f_sky
        self.masks_f_sky.append(self.mask_f_sky)

        # Also want to use the unmasked sky too
        self.masks_f_sky.append(1)

    def read_set_of_masks(self, use_apodization=False):
        """
        Function that reads in the set of masks generated by the generate_set_of_masks function, and stores their
        sky fraction (f_sky) into the class.

        Args:
            use_apodization (bool): Whether to use the apodized version of the mask or not

        Returns:
            None
        """

        for mask_num in range(1, 10):
            try:
                mask = hp.read_map(
                    self.folder_path + 'Masks/Mask' + str(mask_num) + ('_apo' if use_apodization else '') + '.fits')

                f_sky = mask.sum() / mask.size
                f_sky = np.count_nonzero(mask) / mask.size

                print('Fraction of sky allowed through by the dummy mask {msk} after apodizing is: {num:.9f} %'.format(
                    msk=mask_num, num=f_sky * 100))

                self.masks_f_sky.append(f_sky)

            except FileNotFoundError as err:
                # print('Hello', err)
                print(err.filename)
                raise RuntimeError('Oh dear')

        # Append the Euclid mask's f_sky
        self.masks_f_sky.append(self.mask_f_sky)

        # Also want to use the unmasked sky too
        self.masks_f_sky.append(1)

    def custom_mask(self):
        """
        Function to generate a test map that has both the Milky Way and Solar System masked out, but no bright stars yet

        Returns:
            None
        """
        # The N_side parameter for our trial mask
        n_side = 2048
        n_pix = 12 * n_side ** 2

        # Create two blank maps where we initially allow all light through the mask
        map_gal = np.ones(n_pix, dtype=np.bool)
        map_elp = np.ones(n_pix, dtype=np.bool)

        # Create a cylindrical region in galactic coordinates to represent the Milky Way
        region_gal = hp.query_strip(nside=n_side, theta1=1 * np.pi / 3, theta2=2 * np.pi / 3)
        # Create a cylindrical region in ecliptic coordinates to represent the Solar System
        region_elp = hp.query_strip(nside=n_side, theta1=2 * np.pi / 5, theta2=3 * np.pi / 5)

        # Set the regions for the Milky Way and Solar System to zero in the mask
        map_gal[region_gal] = 0
        map_elp[region_elp] = 0

        # Plot both the Milky Way and Solar System regions in their coordinate systems on the same figure
        plt.figure(figsize=(13, 7))

        hp.mollview(map_gal, coord='G', title="Dummy mask of the Milky Way", cmap='seismic', sub=[1, 2, 1],
                    margins=(0.005, 0.0, 0.0, -0.1))
        hp.graticule()

        hp.mollview(map_elp, coord='E', title="Dummy mask of the Solar System", cmap='seismic', sub=[1, 2, 2],
                    margins=(0.0, 0.0, 0.005, -0.1))
        hp.graticule()

        # Combine the two masks into a single mask.
        # Note that we have to rotate the Solar System mask from ecliptic coords to galactic coords
        map_both = np.logical_and(map_gal, hp.rotator.Rotator(coord='EG').rotate_map_pixel(map_elp))

        # Do some manual memory management once got both combined
        del map_gal
        del map_elp

        print('Fraction of sky allowed through by the dummy mask is {num:.2f} %'.format(
            num=100 * map_both.sum() / map_both.size))

        # Plot the resulting mask
        hp.mollview(map_both, coord='G', title="Both masks applied together", cmap='seismic')
        hp.graticule()

        # Also plot the mask using a spherical model
        hp.orthview(map_both, rot=(0, 20, 15), title="Both masks applied together", cmap='seismic')
        hp.graticule()

        ell = np.arange(2, self.ell_max + 1)
        # Now find what the power spectrum of our custom mask and the Euclid mask
        mask_cls = hp.anafast(map_both, lmax=self.ell_max)

        euclid_mask = hp.read_map('./resources/Euclid_masks/Euclid-gal-mask-2048.fits', verbose=False).astype(
                np.bool)

        print('Fraction of sky allowed through by the Euclid mask is {num:.2f} %'.format(
            num=100 * euclid_mask.sum() / euclid_mask.size))

        euclid_mask_cls = hp.anafast(euclid_mask, lmax=self.ell_max)

        # Plot the Cl values for our mask and
        plt.figure(figsize=(13, 7))
        plt.loglog(ell, ell * (ell + 1) * mask_cls[2:] / (2 * np.pi),
                   lw=3, c='tab:blue', label=r'Dummy mask')
        plt.loglog(ell, ell * (ell + 1) * euclid_mask_cls[2:] / (2 * np.pi),
                   lw=3, c='orange', label=r'Euclid mask')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_{\ell} / 2 \pi$')
        plt.title(r'Power spectrum of masks')

        plt.legend()
        plt.tight_layout()

        # * We now want to generate two dummy maps that only have one ell mode activate, with a power of one, and all
        # * other Cl values equal to zero. From here, we then want to apply a mask and recover the power spectrum
        # * seeing how our single Cl mode that we started with bleeds across into other modes.

        # The ell value that we will activate
        ell1 = 50

        # Set up two dummy sets of Cl values that are zero, except for a single ell mode
        dummy_cls1 = np.zeros(self.ell_max + 1)
        dummy_cls1[ell1] = 2 * np.pi / (ell1 * (ell1 + 1))

        # Convert our Cl array into a map
        dummy_map1 = hp.synfast(dummy_cls1, nside=n_side, lmax=self.ell_max)

        # Now mask our our map using the Euclid mask
        dummy_masked_map1 = hp.ma(dummy_map1)
        dummy_masked_map1.mask = np.logical_not(euclid_mask)

        # Mask our map using our existing dummy map
        dummy_masked_map2 = hp.ma(dummy_map1)
        dummy_masked_map2.mask = np.logical_not(map_both)

        # Recover the Cl's from the masked maps
        recov_cls1 = hp.anafast(dummy_masked_map1, lmax=self.ell_max)
        recov_cls2 = hp.anafast(dummy_masked_map2, lmax=self.ell_max)

        # Plot the original and recovered sets of Cl values
        plt.figure(figsize=(11, 6))

        plt.semilogy(ell, ell * (ell + 1) * dummy_cls1[2:] / (2 * np.pi),
                     lw=3, c='tab:blue', label=r'Input value')
        plt.semilogy(ell, ell * (ell + 1) * recov_cls1[2:] / (2 * np.pi),
                     lw=3, c='orange', label=r'Euclid mask')
        plt.semilogy(ell, ell * (ell + 1) * recov_cls2[2:] / (2 * np.pi),
                     lw=3, c='hotpink', label=r'Basic mask', alpha=0.7)

        plt.xlim(left=0, right=100)
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_{\ell} / 2 \pi$')
        plt.title(r'Recovered power spectrum of a single mode with a mask applied')

        plt.legend()
        plt.tight_layout()
        plt.show()

    def theta_phi_only_mask(self):
        """
        Function that investigates masks that only vary in the theta OR the phi direction.

        Returns:
            None
        """

        @functools.lru_cache(maxsize=None)
        def gamma(x):
            """
            Function that implements the Gamma function using mixed-precision floating point calculations and
            function memoization to store the results in memory

            Args:
                x (float): Input value

            Returns:
                mpmath.mp: mpmath's mixed-precision floating point number
            """
            return mpmath.gamma(x)

        def factorial(x):
            """
            Converts the factorial function to a Gamma function call

            Args:
                x (float): Input value

            Returns:
                mpmath.mp: mpmath's mixed-precision floating point number
            """
            return gamma(x + 1)

        def theory_alm_phi(el, m, a_lim):
            """
            Function that returns the squared alm value for a mask that only varies in the phi direction, for a given
            ell, m combination.

            Args:
                el (int): The l value for this alm
                m (int): The m value for this alm
                a_lim (float): Half the width of the mask centred around phi = phi

            Returns:
                float
            """
            alm = (2 + 2 * (-1) ** (el + m)) * mpmath.power(2, 2 * m - 4)
            alm *= (gamma(el / 2) * gamma((el + m + 1) / 2) / (gamma((el + 3) / 2) * factorial((el - m) / 2))) ** 2

            alm *= factorial(el - m) / factorial(el + m)
            alm *= mpmath.sin(m * a_lim) ** 2 / np.pi

            return alm

        def theory_cl_phi(el, a_lim):
            """
            Function that returns the theoretical Cl value for a mask that only depends in the phi direction.

            Args:
                el (int): The l value for this Cl
                a_lim (float): Half the width of the mask centred around phi = pi

            Returns:
                float
            """
            summand = 0

            for m in range(-el, el + 1):
                summand += theory_alm_phi(el, m, a_lim)

            return summand

        n_side = 2048
        n_pix = 12 * n_side ** 2

        # Initiate two masks for our theta and phi cuts which have all pixels blocked initially
        map_theta = np.zeros(n_pix)
        map_phi = np.zeros(n_pix)

        # Region allowed through in a strip centred around theta = pi / 2
        a = np.pi / 3
        b = 2 * np.pi / 3

        region_theta = hp.query_strip(nside=n_side, theta1=a, theta2=b)

        # Set the above pixels to a value of one in the map
        map_theta[region_theta] = 1

        # Convert our theta mask into a set of Cl values
        mask_theta_cls = hp.anafast(map_theta, lmax=self.ell_max)

        ells = np.arange(2, self.ell_max + 1)

        # Use our analytic expressions to get approximate values of the theoretically predicted Cl values
        theory_cls = np.pi * (-scispec.lpmv(0, ells - 1, np.cos(a)) + scispec.lpmv(0, ells + 1, np.cos(a))
                              + scispec.lpmv(0, ells - 1, np.cos(b)) - scispec.lpmv(0, ells + 1, np.cos(b))) ** 2 / (
                             (2 * ells + 1) ** 2)

        # Plot the theory vs numerical Cl values
        plt.figure(figsize=(13, 7))
        plt.loglog(ells, ells * (ells + 1) * theory_cls / (2 * np.pi),
                   lw=3, c='orange', label=r'Theory vals')

        plt.loglog(ells, ells * (ells + 1) * mask_theta_cls[2:] / (2 * np.pi),
                   lw=2, c='tab:blue', label=r'Recov from mask')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_{\ell} / 2 \pi$')
        plt.title(r'Power spectrum of the simple $\vartheta$ mask')

        plt.legend()
        plt.tight_layout()

        # Now make the same plot, but only for the even ells
        plt.figure(figsize=(13, 7))

        ells_even = ells[::2]
        plt.loglog(ells_even, ells_even * (ells_even + 1) * theory_cls[::2] / (2 * np.pi),
                   lw=3, c='orange', label=r'Theory vals')

        plt.loglog(ells_even, ells_even * (ells_even + 1) * mask_theta_cls[2::2] / (2 * np.pi),
                   lw=1, c='tab:blue', label=r'Recov from mask')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_{\ell} / 2 \pi$')
        plt.title(r'Power spectrum of the simple $\vartheta$ mask for even $\ell$ only')

        plt.legend()
        plt.tight_layout()

        # Now plot the ratio of our even Cl's
        plt.figure(figsize=(13, 7))
        plt.semilogx(ells_even, mask_theta_cls[2::2] / theory_cls[::2],
                     lw=3, c='tab:blue', label=r'Ratio')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$C_{\ell}^{\textrm{Num}} / C_{\ell}^{\textrm{Th}} $')
        plt.title(r'Ratio of the numerical to theoretical values for the simple $\vartheta$ mask')

        plt.legend()
        plt.tight_layout()

        plt.show()

        # Create our phi mask, where we only allow data through between phi = pi - a and phi = pi + a
        a = np.pi / 3
        
        for idx in range(len(map_phi)):
            theta, phi = hp.pix2ang(2048, idx)
            if (phi > np.pi - a) and (phi < np.pi + a):
                map_phi[idx] = 1

        # View the mask, but with phi=pi in the centre of the map
        hp.mollview(map_phi, rot=[0, 180, 0], cmap='seismic', title=r'Simple mask in the $\phi$ plane')
        hp.graticule(c='white')
        plt.show()

        # Convert our mask to a set of numerical Cl's
        phi_cls = hp.anafast(map_phi, lmax=self.ell_max)

        # Use the theoretical functions to get a set of theory Cl's
        theory_vals = [theory_cl_phi(i, a) for i in range(2, 2000 + 1)]

        # Plot the numerical vs theory Cl's
        plt.figure(figsize=(13, 7))

        plt.loglog(ells, ells * (ells + 1) * phi_cls[2:] / (2 * np.pi),
                   lw=3, c='tab:blue', label=r'Recov from mask')

        plt.loglog(ells, ells * (ells + 1) * theory_vals / (2 * np.pi),
                   lw=1.5, c='Orange', label=r'Theory vals')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_{\ell} / 2 \pi$')
        plt.title(r'Power spectrum of the simple $\phi$ mask')

        plt.legend()
        plt.tight_layout()

        # Plot the ratio
        plt.figure(figsize=(13, 7))

        plt.semilogx(ells, phi_cls[2:] / theory_vals, lw=2, c='tab:blue', label=r'Ratio')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$C_{\ell}^{\textrm{Num}} / C_{\ell}^{\textrm{Th}} $')
        plt.title(r'Ratio of the numerical to theoretical values for the simple $\phi$ mask')

        plt.legend()
        plt.tight_layout()

        plt.show()

    def noise_simulations(self):
        """
        Function that generates random Gaussian noise due to the intrinsic ellipticities of galaxies in each pixel.

        Returns:
            None
        """

        # Read in the convergence map that has been previously calculated
        converg_map = hp.read_map(self.folder_path + 'Output-poisson-map-f2z2.fits', verbose=False, field=None)

        # Convert this map to a set of Cl values
        converg_cls = hp.anafast(converg_map, lmax=self.ell_max)

        # Show the map
        hp.mollview(converg_map, cmap=make_planck_colour_map(), title='Converg data')
        hp.graticule(verbose=False, alpha=0.8)

        # The N_side parameter of the map generated above
        n_side = 2048

        # Convert our N_side parameter to the number of pixels in the map
        n_pix = 12 * n_side ** 2

        intrinsic_gal_ellip = 0.21  # The standard deviation of the intrinsic galaxy ellipticity distribution

        avg_gal_den = 30  # This is the average surface galaxy density in [num gals / arc min^2]
        area_per_pix = 1.49E8 / n_pix  # This is the total area in arcmin^2 divided by the number of pixels
        num_gal_per_pix = avg_gal_den * area_per_pix
        print('The average number of galaxies per pixel is {num:.3f}'.format(num=num_gal_per_pix))
        print('The standard deviation of our shape noise per pixel is {num:.3f}'
              .format(num=intrinsic_gal_ellip / np.sqrt(num_gal_per_pix)))

        # Generate a map which is just Gaussian random noise, with a mean of zero and st.dev. appropriate
        random_map1 = np.random.normal(loc=0, scale=intrinsic_gal_ellip / np.sqrt(num_gal_per_pix), size=n_pix)

        # Compute the Cl values for our noise map
        cl_noise = hp.anafast(random_map1, lmax=self.ell_max)

        # Compute what the expected Cl VALUE (singular) is for the shape noise.
        # Here, it is the intrinsic ellipticity variance divided by the number density of galaxies (in radians^2)
        theory_cl_noise = intrinsic_gal_ellip ** 2 / (avg_gal_den / (sciconst.arcminute ** 2))

        # Show our random map
        hp.mollview(random_map1, title="Gasussian random noise in each pixel", cmap=make_planck_colour_map())
        hp.graticule(verbose=False, alpha=0.8)

        # Add the noise to our existing convergence map
        converg_w_noise_map = converg_map + random_map1

        # Compute the Cl values for our map with noise
        converg_cls_w_noise = hp.anafast(converg_w_noise_map, lmax=self.ell_max)
        ell = np.arange(2, self.ell_max + 1)

        # Plot the various Cl's
        plt.figure(figsize=(13, 7))
        plt.loglog(ell, ell * (ell + 1) * converg_cls[2:] / (2 * np.pi),
                   lw=3, c='tab:blue', label=r'$C_\ell$ recovered from map')
        plt.loglog(ell, ell * (ell + 1) * cl_noise[2:] / (2 * np.pi),
                   lw=3, c='orange', label=r'Random shape noise')
        plt.loglog(ell, ell * (ell + 1) * converg_cls_w_noise[2:] / (2 * np.pi),
                   lw=2, c='tab:green', label=r'Signal + rand shp noise')
        plt.loglog(ell, ell * (ell + 1) * self.raw_c_ells['W4xW4'][2:] / (2 * np.pi) * np.sqrt(2 / (2 * ell + 1)),
                   lw=2, c='tab:pink', label=r'Cosmic variance')
        plt.loglog(ell, ell * (ell + 1) * self.raw_c_ells['W4xW4'][2:] / (2 * np.pi),
                   lw=2, c='tab:cyan', label=r'Input $C_\ell$')
        plt.loglog(ell, ell * (ell + 1) * theory_cl_noise / (2 * np.pi),
                   lw=2, c='gold', label=r'Theoretical noise')
        plt.loglog(ell, ell * (ell + 1) * (self.raw_c_ells['W4xW4'][2:] + theory_cl_noise) / (2 * np.pi),
                   lw=1.25, c='lime', label=r'Theory $C_\ell$ + noise')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell(\ell+1)C_{\ell} / 2 \pi$')
        plt.title('Power spectrum of lensing with and without Gaussian random noise')

        plt.legend()
        plt.tight_layout()
        plt.show()

        import sys
        sys.exit()

    @staticmethod
    def manually_making_map():
        ell = 75

        alms = np.zeros(ell + 1, dtype=np.complex64)

        rng = np.random.default_rng()

        alms[0] = rng.standard_normal()

        for m in range(1, ell + 1):
            alms[m] = (rng.standard_normal() + 1j * rng.standard_normal())

        n_side = 256
        n_pix = 12 * n_side ** 2

        map1 = np.zeros(n_pix)

        for idx in range(len(map1)):
            theta, phi = hp.pix2ang(n_side, idx)
            for m in range(0, ell + 1):
                if m == 0:
                    map1[idx] += alms[m] * scispec.sph_harm(m, ell, phi, theta)
                else:
                    map1[idx] += 2 * np.real(alms[m] * scispec.sph_harm(m, ell, phi, theta))

        hp.mollview(map1, cmap='seismic')
        plt.show()

    @staticmethod
    def global_minima(x_vals, y_vals, return_y=False):
        """
        Function that computes the global minima of some data, that is assumed to follow y = f(x).

        Args:
            x_vals (ndarray): Array of x values
            y_vals (ndarray): Array of y values evaluated at the x_vals
            return_y (bool): Optional bool if the y value at the global minima should be returned too

        Returns:
            Either list(double, double) or double depending on return_y
        """

        # Wrap in a try-except as sometimes Numpy returns an error finding the minima
        try:
            # First, interpolate the provided x and y values
            # Note that we use a spline of order 4, which guarantees that the derivative is a cubic spline
            spline = sciinterp.InterpolatedUnivariateSpline(x_vals, y_vals, k=4)

            # Now find where the derivative of the spline is equal to zero, i.e. a local minima or maxima
            crit_pts = spline.derivative().roots()

            # Now evaluate the spline at the local minima/maxima to find which index is the true global minima
            min_index = np.argmin(spline(crit_pts))

            # Find what the x value is at this global minima
            crit_x = crit_pts[min_index]

            # Also find the y value
            crit_y = spline(crit_x)

            # Return either the x and y values at this point, or just the x value
            return [crit_x, crit_y] if return_y else crit_x

        except ValueError:
            print('NumPy cannot find zeros of spline, returning zero')
            return 0

    def simple_likelihood_as(self):
        """
        Function that implements a *very* simple likelihood calculation that allows us to estimate the value of A_s
        from a convergence map

        Returns:
            None
        """
        lmax = 2000
        ells = np.arange(2, lmax + 1)

        # The N_side parameter of the generated maps & masks
        n_side = 2048

        # Convert our N_side parameter to the number of pixels in the map
        n_pix = 12 * n_side ** 2

        intrinsic_gal_ellip = 0.21  # The standard deviation of the intrinsic galaxy ellipticity distribution

        avg_gal_den = 30  # This is the average surface galaxy density in [num gals / arc min^2]
        area_per_pix = 1.49E8 / n_pix  # This is the total area in arcmin^2 divided by the number of pixels
        num_gal_per_pix = avg_gal_den * area_per_pix

        # Generate random Gaussian noise that will be added to our maps
        random_noise = np.random.normal(loc=0, scale=intrinsic_gal_ellip / np.sqrt(num_gal_per_pix), size=n_pix)

        # Compute what the expected Cl VALUE (singular) is for the shape noise.
        theory_cl_noise = intrinsic_gal_ellip ** 2 / (avg_gal_den / (sciconst.arcminute ** 2))

        theory_cl_noise = ells * (ells + 1) * theory_cl_noise / (2 * np.pi)

        # ! theory_cl_noise = np.zeros(lmax-1)

        # Initiate a LCDM cosmology
        params = camb.CAMBparams()
        params.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0.12)
        params.InitPower.set_params(As=2.1E-9, ns=0.96)
        params.set_for_lmax(lmax, lens_potential_accuracy=1)

        params.Want_CMB = False
        params.WantCls = True
        params.NonLinear = camb.model.NonLinear_both

        # We want to evaluate the lensing power spectrum at z=2 only, for the moment.
        params.SourceWindows = [GaussianSourceWindow(redshift=2, source_type='lensing', sigma=0.05)]

        # The range of A_s values that we want to compute the power spectrum at
        a_s_center = 2.1615E-9
        a_s_vals = np.linspace(a_s_center - 5E-12, a_s_center + 5E-12, 10)

        # List which our Cl values will get stored into
        cl_vals = []

        print('Evaluating power spectrum on A_s values')

        # Go through each A_s value, compute the lensing power spec, and save to the list.
        for a_s in a_s_vals:
            params.InitPower.set_params(As=a_s, ns=0.96)

            results = camb.get_results(params)

            cls = results.get_source_cls_dict()

            cl_vals.append(cls['W1xW1'][2:lmax + 1])

        # Now we want to spline the power spectrum at each ell to get this as a function of A_s.
        # This allows us to predict the power spectra at arbitrary ell.
        splines = []

        print('Evaluating splines')

        for ell_idx in range(len(ells)):
            splines.append(sciinterp.InterpolatedUnivariateSpline(a_s_vals,
                                                                  [cl_vals[a_s_idx][ell_idx] for a_s_idx in
                                                                   range(len(a_s_vals))]))

        print('Reading in convergence & shear maps')

        # Read in the unmasked convergence map that was created with A_s = 2.25E-9 and m_nu = 0.12
        converg_map, gamma1_map, gamma2_map = hp.read_map(self.folder_path + 'KappaGammaMap-f2z2.fits', verbose=True,
                                                          field=None)

        # Add the random shape noise to our maps
        converg_map_noise = converg_map + random_noise
        gamma1_map_noise = gamma1_map + random_noise
        gamma2_map_noise = gamma2_map + random_noise

        # Convert the unmasked maps to C_ells
        converg_cls, shear_cls = hp.anafast([converg_map, gamma1_map, gamma2_map], lmax=lmax, nspec=2)[:, 2:]
        converg_cls = ells * (ells + 1) * converg_cls / (2 * np.pi)
        shear_cls = ells * (ells + 1) * shear_cls / (2 * np.pi)

        converg_cls_noise, shear_cls_noise = hp.anafast([converg_map_noise, gamma1_map_noise, gamma2_map_noise],
                                                        lmax=lmax, nspec=2)[:, 2:]
        converg_cls_noise = ells * (ells + 1) * converg_cls_noise / (2 * np.pi)
        shear_cls_noise = ells * (ells + 1) * shear_cls_noise / (2 * np.pi)

        print('Reading in the Euclid mask')
        # Read in the Euclid-like mask
        euclid_mask = hp.read_map('./resources/Euclid_masks/Euclid-gal-mask-2048.fits', verbose=False).astype(
            np.bool)
        sky_fraction = euclid_mask.sum() / euclid_mask.size

        print('Masking the convergence & shear maps')

        # Mask the convergence map
        masked_converg_map = hp.ma(converg_map)
        masked_converg_map.mask = np.logical_not(euclid_mask)

        # Mask the shear maps
        masked_gamma1_map = hp.ma(gamma1_map)
        masked_gamma2_map = hp.ma(gamma2_map)
        masked_gamma1_map.mask = np.logical_not(euclid_mask)
        masked_gamma2_map.mask = np.logical_not(euclid_mask)

        # Mask the convergence map with noise
        masked_converg_map_noise = hp.ma(converg_map_noise)
        masked_converg_map_noise.mask = np.logical_not(euclid_mask)

        # Mask the shear maps with noise
        masked_gamma1_map_noise = hp.ma(gamma1_map_noise)
        masked_gamma2_map_noise = hp.ma(gamma2_map_noise)
        masked_gamma1_map_noise.mask = np.logical_not(euclid_mask)
        masked_gamma2_map_noise.mask = np.logical_not(euclid_mask)

        print('Turning maps into power spectra')

        # Obtain masked Cl's for both masks
        masked_converg_cls, masked_shear_cls = hp.anafast([masked_converg_map, masked_gamma1_map, masked_gamma2_map],
                                                          lmax=lmax, nspec=2)[:, 2:]
        masked_converg_cls = ells * (ells + 1) * masked_converg_cls / (2 * np.pi) / sky_fraction
        masked_shear_cls = ells * (ells + 1) * masked_shear_cls / (2 * np.pi) / sky_fraction

        masked_converg_cls_noise, masked_shear_cls_noise = hp.anafast(
            [masked_converg_map_noise, masked_gamma1_map_noise, masked_gamma2_map_noise], lmax=lmax, nspec=2)[:, 2:]
        masked_converg_cls_noise = ells * (ells + 1) * masked_converg_cls_noise / (2 * np.pi) / sky_fraction
        masked_shear_cls_noise = ells * (ells + 1) * masked_shear_cls_noise / (2 * np.pi) / sky_fraction

        print('Plotting power spectra')

        # Compute the "true" values for the convergence and shear signals, with no noise and masking
        true_converg = np.array([splines[ell](2.1E-9) for ell in range(len(ells))])
        true_shear = (ells + 2) * (ells - 1) / (ells * (ells + 1)) * true_converg

        # Plot the different convergence signals
        plt.figure(figsize=(11, 6))
        plt.loglog(ells, true_converg,
                   lw=1.5, ls='-', c='cornflowerblue', label=r'True convergence')
        plt.loglog(ells, converg_cls,
                   lw=1, ls='-', c='tab:blue', label=r'No mask no noise')
        plt.loglog(ells, converg_cls_noise,
                   lw=1, ls='-', c='orange', label=r'No mask with noise')
        plt.loglog(ells, masked_converg_cls,
                   lw=1, ls='-', c='hotpink', label=r'With mask no noise')
        plt.loglog(ells, masked_converg_cls_noise,
                   lw=1, ls='-', c='purple', label=r'With mask and noise')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell^{\kappa \kappa} / 2 \pi$')
        plt.legend()
        plt.tight_layout()

        # Plot the shear signals
        plt.figure(figsize=(11, 6))
        plt.loglog(ells, true_shear,
                   lw=1.5, ls='-', c='cornflowerblue', label=r'True shear')
        plt.loglog(ells, shear_cls,
                   lw=1, ls='-', c='tab:blue', label=r'No mask no noise')
        plt.loglog(ells, shear_cls_noise,
                   lw=1, ls='-', c='orange', label=r'No mask with noise')
        plt.loglog(ells, masked_shear_cls,
                   lw=1, ls='-', c='hotpink', label=r'With mask no noise')
        plt.loglog(ells, masked_shear_cls_noise,
                   lw=1, ls='-', c='purple', label=r'With mask and noise')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell^{EE} / 2 \pi$')
        plt.legend()
        plt.tight_layout()

        # Plot the true values of the convergence and shear
        plt.figure(figsize=(11, 6))
        plt.loglog(ells, true_converg,
                   lw=1, ls='-', c='tab:blue', label=r'True convergence')
        plt.loglog(ells, true_shear,
                   lw=1, ls='-', c='orange', label=r'True shear')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell^{ii} / 2 \pi$')
        plt.legend()
        plt.tight_layout()

        plt.show()

        # Now we want to evaluate the likelihood for a range of A_s values to find the maximum-likelihood value
        a_s_vals = np.linspace(a_s_center - 3.5E-12, a_s_center + 2.5E-12, 250)

        # Likelihood for our unmasked, Euclid mask, and custom mask
        likelihoods1_converg = []
        likelihoods2_converg = []
        likelihoods3_converg = []
        likelihoods4_converg = []

        likelihoods1_shear = []
        likelihoods2_shear = []
        likelihoods3_shear = []
        likelihoods4_shear = []

        for a_s in a_s_vals:
            cl_theory = np.array([splines[ell](a_s) for ell in range(len(ells))])

            # * Compute what the log-likelihood is for the four cases:
            # Unmasked map with no noise
            log_lik = np.sum((2 * ells + 1) * (np.log(cl_theory) + converg_cls / cl_theory))

            # Unmasked map with noise
            log_lik2 = np.sum((2 * ells + 1) * (np.log(cl_theory + theory_cl_noise) +
                                                converg_cls_noise / (cl_theory + theory_cl_noise)))

            # Masked map without noise
            log_lik3 = np.sum((2 * ells + 1) * (np.log(cl_theory) + masked_converg_cls / cl_theory))

            # Masked map with noise
            log_lik4 = np.sum((2 * ells + 1) * (np.log(cl_theory + theory_cl_noise)
                                                + masked_converg_cls_noise / (cl_theory + theory_cl_noise)))

            likelihoods1_converg.append(log_lik)
            likelihoods2_converg.append(log_lik2)
            likelihoods3_converg.append(log_lik3)
            likelihoods4_converg.append(log_lik4)

            # Now turn the convergence Cl's into EE
            cl_theory = (ells + 2) * (ells - 1) / (ells * (ells + 1)) * cl_theory

            log_lik = np.sum((2 * ells + 1) * (np.log(cl_theory) + converg_cls / cl_theory))

            # Unmasked map with noise
            log_lik2 = np.sum((2 * ells + 1) * (np.log(cl_theory + theory_cl_noise) +
                                                converg_cls_noise / (cl_theory + theory_cl_noise)))

            # Masked map without noise
            log_lik3 = np.sum((2 * ells + 1) * (np.log(cl_theory) + masked_converg_cls / cl_theory))

            # Masked map with noise
            log_lik4 = np.sum((2 * ells + 1) * (np.log(cl_theory + theory_cl_noise)
                                                + masked_converg_cls_noise / (cl_theory + theory_cl_noise)))

            likelihoods1_shear.append(log_lik)
            likelihoods2_shear.append(log_lik2)
            likelihoods3_shear.append(log_lik3)
            likelihoods4_shear.append(log_lik4)

        # Now compute the global minima of the likelihoods1_converg, to find the maximum-likelihood value of A_s
        cr_pts1_cn = self.global_minima(a_s_vals, likelihoods1_converg)
        cr_pts2_cn = self.global_minima(a_s_vals, likelihoods2_converg)
        cr_pts3_cn = self.global_minima(a_s_vals, likelihoods3_converg)
        cr_pts4_cn = self.global_minima(a_s_vals, likelihoods4_converg)

        print(f'Unmasked no noise convergence \t A_s: {cr_pts1_cn:.5e}')
        print(f'Unmasked with noise convergence  A_s: {cr_pts2_cn:.5e}')
        print(f'Masked no noise convergence \t A_s: {cr_pts3_cn:.5e}')
        print(f'Masked with noise convergence \t A_s: {cr_pts4_cn:.5e}\n')

        cr_pts1_sh = self.global_minima(a_s_vals, likelihoods1_shear)
        cr_pts2_sh = self.global_minima(a_s_vals, likelihoods2_shear)
        cr_pts3_sh = self.global_minima(a_s_vals, likelihoods3_shear)
        cr_pts4_sh = self.global_minima(a_s_vals, likelihoods4_shear)

        print(f'Unmasked no noise shear \t A_s: {cr_pts1_sh:.5e}')
        print(f'Unmasked with noise shear \t A_s: {cr_pts2_sh:.5e}')
        print(f'Masked no noise shear \t\t A_s: {cr_pts3_sh:.5e}')
        print(f'Masked with noise shear \t A_s: {cr_pts4_sh:.5e}\n')

        # * Plot the likelihood as a function of A_s
        # Note that we have to used a zoomed in sub-plot as the lines aren't separated by much

        """
        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        from mpl_toolkits.axes_grid1.inset_locator import mark_inset

        fig, ax = plt.subplots(figsize=(13, 7))
        ax.scatter(a_s_vals, likelihoods, c='tab:blue', label='Unmasked')
        ax.scatter(a_s_vals, likelihoods_2, c='orange', label='Masked')
        ax.scatter(a_s_vals, likelihoods_3, c='hotpink', label='Masked2')

        ax2 = zoomed_inset_axes(ax, 7, loc=1, borderpad=1)

        ax2.scatter(a_s_vals, likelihoods, c='tab:blue', label='Unmasked')
        ax2.scatter(a_s_vals, likelihoods_2, c='orange', label='Euclid mask')
        ax2.scatter(a_s_vals, likelihoods_3, c='hotpink', label='Mask 2')

        ax2.set_xlim(2.05e-9, 2.45e-9)
        ax2.set_ylim(-3.105e7, -3.08e7)

        mark_inset(ax, ax2, loc1=2, loc2=4, fc="none", ec="0.25")

        ax.axvline(x=cr_pts, c='tab:blue')
        ax.axvline(x=cr_pts2, c='orange')
        ax.axvline(x=cr_pts3, c='hotpink')
        ax.axvline(x=2.25E-9, c='tab:cyan', ls='--')

        ax2.axvline(x=cr_pts, c='tab:blue')
        ax2.axvline(x=cr_pts2, c='orange')
        ax2.axvline(x=cr_pts3, c='hotpink')
        ax2.axvline(x=2.25E-9, c='tab:cyan', ls='--')

        ax.set_title(r'Log-likelihood as a function of $A_\textsc{s}$ for masked and unmasked maps')
        ax.set_xlabel(r'$A_\textsc{s}$')
        ax.set_ylabel(r'$\propto -\ln \mathcal{L}$')

        ax.legend(loc='lower right')
        fig.tight_layout()
        """

        likelihoods1_converg = np.array(likelihoods1_converg)
        likelihoods2_converg = np.array(likelihoods2_converg)
        likelihoods3_converg = np.array(likelihoods3_converg)
        likelihoods4_converg = np.array(likelihoods4_converg)

        likelihoods1_shear = np.array(likelihoods1_shear)
        likelihoods2_shear = np.array(likelihoods2_shear)
        likelihoods3_shear = np.array(likelihoods3_shear)
        likelihoods4_shear = np.array(likelihoods4_shear)

        min_val1_cn = np.min(likelihoods1_converg)
        a_s_center1_cn = a_s_vals[np.argmin(likelihoods1_converg)]

        min_val2_cn = np.min(likelihoods2_converg)
        a_s_center2_cn = a_s_vals[np.argmin(likelihoods2_converg)]

        min_val3_cn = np.min(likelihoods3_converg)
        a_s_center3_cn = a_s_vals[np.argmin(likelihoods3_converg)]

        min_val4_cn = np.min(likelihoods4_converg)
        a_s_center4_cn = a_s_vals[np.argmin(likelihoods4_converg)]

        # Now do the same for shear
        min_val1_sh = np.min(likelihoods1_shear)
        a_s_center1_sh = a_s_vals[np.argmin(likelihoods1_shear)]

        min_val2_sh = np.min(likelihoods2_shear)
        a_s_center2_sh = a_s_vals[np.argmin(likelihoods2_shear)]

        min_val3_sh = np.min(likelihoods3_shear)
        a_s_center3_sh = a_s_vals[np.argmin(likelihoods3_shear)]

        min_val4_sh = np.min(likelihoods4_shear)
        a_s_center4_sh = a_s_vals[np.argmin(likelihoods4_shear)]

        # First plot the log-likelihood values
        fig, ax = plt.subplots(figsize=(13, 7))

        ax.plot(a_s_vals, likelihoods1_converg - min_val1_cn, c='tab:blue', label='Unmasked no noise')
        ax.plot(a_s_vals, likelihoods2_converg - min_val2_cn, c='orange', label='Unmasked with noise')
        ax.plot(a_s_vals, likelihoods3_converg - min_val3_cn, c='hotpink', label='Masked no noise', ls='--')
        ax.plot(a_s_vals, likelihoods4_converg - min_val4_cn, c='purple', label='Masked with noise', ls='--')

        ax.set_xlabel(r'$A_\textsc{s}$')
        ax.set_ylabel(r'$\ln \mathcal{L} - \textrm{Min}[\ln \mathcal{L}]$')
        ax.set_title(r'Normalised Log-likelihoods values for convergence-$\kappa\kappa$')
        plt.legend()
        plt.tight_layout()

        # Plot the log-likelihood values for shear
        fig, ax = plt.subplots(figsize=(13, 7))

        ax.plot(a_s_vals, likelihoods1_shear - min_val1_sh, c='tab:blue', label='Unmasked no noise')
        ax.plot(a_s_vals, likelihoods2_shear - min_val2_sh, c='orange', label='Unmasked with noise')
        ax.plot(a_s_vals, likelihoods3_shear - min_val3_sh, c='hotpink', label='Masked no noise', ls='--')
        ax.plot(a_s_vals, likelihoods4_shear - min_val4_sh, c='purple', label='Masked with noise', ls='--')

        ax.set_xlabel(r'$A_\textsc{s}$')
        ax.set_ylabel(r'$\ln \mathcal{L} - \textrm{Min}[\ln \mathcal{L}]$')
        ax.set_title(r'Normalised Log-likelihoods values for shear-$EE$')
        plt.legend()
        plt.tight_layout()

        # Now plot the normalised likelihood for convergence
        fig, ax = plt.subplots(figsize=(13, 7))

        ax.plot(a_s_vals - a_s_center1_cn, np.exp(likelihoods1_converg - min_val1_cn), c='tab:blue',
                label='Unmasked no noise')
        ax.plot(a_s_vals - a_s_center2_cn, np.exp(likelihoods2_converg - min_val2_cn), c='orange',
                label='Unmasked with noise')
        ax.plot(a_s_vals - a_s_center3_cn, np.exp(likelihoods3_converg - min_val3_cn), c='hotpink',
                label='Masked no noise', ls='--')
        ax.plot(a_s_vals - a_s_center4_cn, np.exp(likelihoods4_converg - min_val4_cn), c='purple',
                label='Masked with noise', ls='--')

        ax.set_xlim(left=-2E-12, right=2E-12)
        ax.set_ylim(bottom=0.8, top=9.1)

        ax.set_xlabel(r'$A_\textsc{s} - A_\textsc{s}^\textsc{ml}$')
        ax.set_ylabel(r'$\sim -\mathcal{L}$')
        ax.set_title(r'Normalised likelihoods where $\textrm{Min}[\mathcal{L}] = 1$ for convergence-$\kappa\kappa$')
        plt.legend()
        plt.tight_layout()

        # Now plot the normalised likelihood for shear
        fig, ax = plt.subplots(figsize=(13, 7))

        ax.plot(a_s_vals - a_s_center1_sh, np.exp(likelihoods1_shear - min_val1_sh), c='tab:blue',
                label='Unmasked no noise')
        ax.plot(a_s_vals - a_s_center2_sh, np.exp(likelihoods2_shear - min_val2_sh), c='orange',
                label='Unmasked with noise')
        ax.plot(a_s_vals - a_s_center3_sh, np.exp(likelihoods3_shear - min_val3_sh), c='hotpink',
                label='Masked no noise', ls='--')
        ax.plot(a_s_vals - a_s_center4_sh, np.exp(likelihoods4_shear - min_val4_sh), c='purple',
                label='Masked with noise', ls='--')

        ax.set_xlim(left=-2E-12, right=2E-12)
        ax.set_ylim(bottom=0.8, top=9.1)

        ax.set_xlabel(r'$A_\textsc{s} - A_\textsc{s}^\textsc{ml}$')
        ax.set_ylabel(r'$\sim -\mathcal{L}$')
        ax.set_title(r'Normalised likelihoods where $\textrm{Min}[\mathcal{L}] = 1$ for shear-$EE$')
        plt.legend()
        plt.tight_layout()

        plt.show()

    def simple_likelihood_mnu(self):
        """
        Function that implements a *very* simple likelihood calculation that allows us to estimate the value of m_nu
        from a convergence map

        Returns:
            None
        """
        lmax = 2000
        ells = np.arange(2, lmax + 1)

        # Initiate a LCDM cosmology
        params = camb.CAMBparams()
        params.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0.12)
        params.InitPower.set_params(As=2.25E-9, ns=0.96)
        params.set_for_lmax(lmax, lens_potential_accuracy=1)

        params.Want_CMB = False
        params.WantCls = True
        params.NonLinear = camb.model.NonLinear_both

        # We want to evaluate the lensing power spectrum at z=2 only, for the moment.
        params.SourceWindows = [GaussianSourceWindow(redshift=2, source_type='lensing', sigma=0.05)]

        # The range of m_nu values that we want to compute the power spectrum at
        m_nu_vals = np.linspace(0, 0.25, 10)

        # List which our Cl values will get stored into
        cl_vals = []

        # Go through each m_nu value, compute the lensing power spec, and save to the list.
        for m_nu in m_nu_vals:
            params.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=m_nu)

            results = camb.get_results(params)

            cls = results.get_source_cls_dict()

            cl_vals.append(cls['W1xW1'][2:lmax + 1])

        # Now we want to spline the power spectrum at each ell to get this as a function of m_nu.
        # This allows us to predict the power spectra at arbitrary ell.
        splines = []

        for ell_idx in range(len(ells)):
            splines.append(sciinterp.InterpolatedUnivariateSpline(m_nu_vals,
                                                                  [cl_vals[m_nu_idx][ell_idx] for m_nu_idx in
                                                                   range(len(m_nu_vals))]))

        # Read in the unmasked convergence map that was created with A_s = 2.25E-9 and m_nu = 0.12
        converg_map = hp.read_map(self.folder_path + 'Output-poisson-map-2-f2z2.fits', verbose=False, field=None)

        # Convert to C_ells
        converg_cls = np.array(hp.anafast(converg_map, lmax=lmax)[2:])
        converg_cls = ells * (ells + 1) * converg_cls / (2 * np.pi)

        # Read in the Euclid-like mask
        euclid_mask = hp.read_map('./resources/Euclid_masks/Euclid-gal-mask-2048.fits', verbose=False).astype(
            np.bool)
        sky_fraction = euclid_mask.sum() / euclid_mask.size

        # Read in our mask that has a crazy small f_sky
        mask2 = hp.read_map(self.folder_path + 'Masks/Mask1_2048.fits', verbose=False).astype(np.bool)
        sky_fraction2 = mask2.sum() / mask2.size

        masked_map = hp.ma(converg_map)
        masked_map.mask = np.logical_not(euclid_mask)

        masked_map2 = hp.ma(converg_map)
        masked_map2.mask = np.logical_not(mask2)

        # Obtain masked Cl's for both masks
        masked_cls = np.array(hp.anafast(masked_map, lmax=lmax)[2:])
        masked_cls = ells * (ells + 1) * masked_cls / (2 * np.pi) / sky_fraction

        masked_cls2 = np.array(hp.anafast(masked_map2, lmax=lmax)[2:])
        masked_cls2 = ells * (ells + 1) * masked_cls2 / (2 * np.pi) / sky_fraction2

        plt.loglog(ells, converg_cls,
                   lw=1, ls='-', c='tab:blue', label=r'$C_\ell$ recovered from full map')
        plt.loglog(ells, masked_cls,
                   lw=1, ls='-', c='orange', label=r'$C_\ell$ with Euclid mask')
        plt.loglog(ells, masked_cls2,
                   lw=1, ls='-', c='hotpink', label=r'$C_\ell$ with mask 2')
        plt.show()

        # Now we want to evaluate the likelihood for a range of A_s values to find the maximum-likelihood value
        m_nu_vals = np.linspace(0, 0.25, 100)

        # Likelihood for our unmasked, Euclid mask, and custom mask
        likelihoods = []
        likelihoods_2 = []
        likelihoods_3 = []

        for m_nu in m_nu_vals:
            cl_theory = np.array([splines[ell](m_nu) for ell in range(len(ells))])

            # Compute what the log-likelihood is for the three cases
            log_lik = np.sum((2 * ells + 1) * (np.log(cl_theory) + converg_cls / cl_theory))
            log_lik2 = np.sum((2 * ells + 1) * (np.log(cl_theory) + masked_cls / cl_theory))
            log_lik3 = np.sum((2 * ells + 1) * (np.log(cl_theory) + masked_cls2 / cl_theory))

            likelihoods.append(log_lik)
            likelihoods_2.append(log_lik2)
            likelihoods_3.append(log_lik3)

        # Now compute the global minima of the likelihoods, to find the maximum-likelihood value of m_nu
        cr_pts = self.global_minima(m_nu_vals, likelihoods)
        cr_pts2 = self.global_minima(m_nu_vals, likelihoods_2)
        cr_pts3 = self.global_minima(m_nu_vals, likelihoods_3)

        print('Unmasked m_nu: ', cr_pts)
        print('Masked m_nu: ', cr_pts2)
        print('Masked2 m_nu: ', cr_pts3)

        # * Plot the likelihood as a function of m_nu
        # Note that we don't have to have the zoomed in sub-plot here, as separation isn't as much

        from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
        from mpl_toolkits.axes_grid1.inset_locator import mark_inset

        fig, ax = plt.subplots(figsize=(13, 7))
        ax.scatter(m_nu_vals, likelihoods, c='tab:blue', label='Unmasked')
        ax.scatter(m_nu_vals, likelihoods_2, c='orange', label='Masked')
        ax.scatter(m_nu_vals, likelihoods_3, c='hotpink', label='Masked2')

        ax.axvline(x=cr_pts, c='tab:blue')
        ax.axvline(x=cr_pts2, c='orange')
        ax.axvline(x=cr_pts3, c='hotpink')
        ax.axvline(x=0.12, c='tab:cyan', ls='--')

        """
        ax2 = zoomed_inset_axes(ax, 7, loc=1, borderpad=1)

        ax2.scatter(m_nu_vals, likelihoods, c='tab:blue', label='Unmasked')
        ax2.scatter(m_nu_vals, likelihoods_2, c='orange', label='Euclid mask')
        ax2.scatter(m_nu_vals, likelihoods_3, c='hotpink', label='Mask 2')

        ax2.set_xlim(2.05e-9, 2.45e-9)
        ax2.set_ylim(-3.105e7, -3.08e7)

        mark_inset(ax, ax2, loc1=2, loc2=4, fc="none", ec="0.25")

        ax2.axvline(x=cr_pts, c='tab:blue')
        ax2.axvline(x=cr_pts2, c='orange')
        ax2.axvline(x=cr_pts3, c='hotpink')
        ax2.axvline(x=0.12, c='tab:cyan', ls='--')
        """

        ax.set_title(r'Log-likelihood as a function of $\Sigma m_{\nu}$ for masked and unmasked maps')
        ax.set_xlabel(r'$\Sigma m_{\nu}$')
        ax.set_ylabel(r'$\propto -\ln \mathcal{L}$')

        ax.legend(loc='lower right')
        fig.tight_layout()
        plt.show()

    def multiple_run_flask_with_shear(self, num_runs):
        """
        Custom function that runs Flask many times, applying ten different masks to the convergence maps to see how
        changing f_sky affects the recovered data

        Args:
            num_runs (int): The number of times that Flask should be run

        Returns:
            None
        """

        if self.masked_cl_out_dir is None:
            raise RuntimeError('In order to save the masked Cl data correctly, the output folder is required! Please'
                               'run the "set_multiple_masked_output" function beforehand.')

        # Time how long the total runs took
        start_time = time.time()

        # Dictionary of which is where we will store the shear values into
        EE_data = {'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                   'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        BB_data = {'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                   'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}

        EB_data = {'Mask1': [], 'Mask2': [], 'Mask3': [], 'Mask4': [], 'Mask5': [], 'Mask6': [],
                   'Mask7': [], 'Mask8': [], 'Mask9': [], 'Mask10': []}
        # Here, 'Mask10' refers to the unmasked values

        # Run Flask the given number of times to generate a power spectra each time
        for run_num in range(num_runs):
            self.run_flask(use_rnd=True)

            for mask_num in range(1, 5):
                cl_df = pd.read_csv(self.folder_path + self.masked_cl_out_dir + '/TEB_Cls_' + str(mask_num) + '.dat',
                                    sep=r'\s+')

                # cl_df.columns = cl_df.columns.str.lstrip()

                # Get the ell values as a numpy array
                ells = cl_df['ell'].to_numpy()

                # Get the c_ell values normalised by ell(ell+1)/2pi, as usual
                # Also normalise by 1 / f_sky for the specific mask in question
                ee_vals = cl_df['EE'].to_numpy() * ells * (ells + 1) / (2 * np.pi) / self.masks_f_sky[
                    mask_num - 1]

                bb_vals = cl_df['BB'].to_numpy() * ells * (ells + 1) / (2 * np.pi) / self.masks_f_sky[
                    mask_num - 1]

                eb_vals = cl_df['EB'].to_numpy() * ells * (ells + 1) / (2 * np.pi) / self.masks_f_sky[
                    mask_num - 1]

                # Create a new dummy dataframe with out C_ell values, with an index determined by the ell values
                ee_df = pd.DataFrame({'Cl': ee_vals}, index=ells)
                bb_df = pd.DataFrame({'Cl': bb_vals}, index=ells)
                eb_df = pd.DataFrame({'Cl': eb_vals}, index=ells)

                EE_data['Mask' + str(mask_num)].append(ee_df.transpose())
                BB_data['Mask' + str(mask_num)].append(bb_df.transpose())
                EB_data['Mask' + str(mask_num)].append(eb_df.transpose())

        # Save results for each mask separately
        for mask_num in range(1, 5):
            EE_data['Mask' + str(mask_num)] = pd.concat(EE_data['Mask' + str(mask_num)], ignore_index=True)
            EE_data['Mask' + str(mask_num)].to_csv(self.folder_path + self.masked_cl_out_dir + '/EE_Cls_Mask' +
                                                   str(mask_num) + '.csv', index=False)

            BB_data['Mask' + str(mask_num)] = pd.concat(BB_data['Mask' + str(mask_num)], ignore_index=True)
            BB_data['Mask' + str(mask_num)].to_csv(self.folder_path + self.masked_cl_out_dir + '/BB_Cls_Mask' +
                                                   str(mask_num) + '.csv', index=False)

            EB_data['Mask' + str(mask_num)] = pd.concat(EB_data['Mask' + str(mask_num)], ignore_index=True)
            EB_data['Mask' + str(mask_num)].to_csv(self.folder_path + self.masked_cl_out_dir + '/EB_Cls_Mask' +
                                                   str(mask_num) + '.csv', index=False)

        # Prints timing statistics
        print('The total time for the runs was {num:.3f} seconds, with an average of {numpersec:.3f} seconds/run'
              .format(num=time.time() - start_time, numpersec=(time.time() - start_time) / num_runs))

        # Finally, go through each MaskedCl output file and delete them once done to tidy up
        for mask_num in range(1, 5):
            try:
                os.remove(self.folder_path + self.masked_cl_out_dir + '/TEB_Cls_' + str(mask_num) + '.dat')
            except OSError as err:
                print("Error: couldn't delete file {file}, error is: {error}".format(file=err.filename,
                                                                                     error=err.strerror))

    def plot_multiple_run_flask_with_shear(self):
        """
        Function to plot the data that was obtained for the shear E/B decomposition using the above function.

        Returns:
            None
        """
        num_maps = 5

        # Temporary variables to construct the dictionaries out of
        tmp_names = ['Mask' + str(num) for num in range(1, num_maps + 1)]
        tmp_value = [np.zeros(self.ell_max - 1) for _ in range(1, num_maps + 1)]
        tmp_dict = dict(zip(tmp_names, tmp_value))

        # Create dictionaries that store the results in
        EE_avg = copy.deepcopy(tmp_dict)
        EE_var = copy.deepcopy(tmp_dict)
        EE_skew = copy.deepcopy(tmp_dict)
        BB_avg = copy.deepcopy(tmp_dict)
        BB_var = copy.deepcopy(tmp_dict)
        BB_skew = copy.deepcopy(tmp_dict)
        EB_rms = copy.deepcopy(tmp_dict)  # Note that for EB we use root-mean-square (RMS) not mean values
        EB_var = copy.deepcopy(tmp_dict)
        EB_skew = copy.deepcopy(tmp_dict)

        # Go through each map to compute statistics
        for mask_num in range(1, num_maps + 1):

            # Print the current working mask number
            print(f'Mask {mask_num}', end='\t', flush=True)

            mask_key = 'Mask' + str(mask_num)

            # * Read in EE data
            data_df = pd.read_csv(
                self.folder_path + self.masked_cl_out_dir + '/EE_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                ell_key = int(label) - 2

                if mask_num == 4:
                    c_ells *= self.mask_f_sky

                # Compute mean
                mean_val = np.mean(c_ells)
                EE_avg[mask_key][ell_key] = mean_val

                # Compute variance
                EE_var[mask_key][ell_key] = np.var(c_ells) / (mean_val * mean_val)

                # Compute skew
                EE_skew[mask_key][ell_key] = scistats.skew(c_ells, bias=True)

            # * Read in BB data
            data_df = pd.read_csv(
                self.folder_path + self.masked_cl_out_dir + '/BB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                ell_key = int(label) - 2

                if mask_num == 4:
                    c_ells *= self.mask_f_sky

                # Compute mean
                mean_val = np.mean(c_ells)
                BB_avg[mask_key][ell_key] = mean_val

                # Compute variance
                BB_var[mask_key][ell_key] = np.var(c_ells) / (mean_val * mean_val)

                # Compute skew
                BB_skew[mask_key][ell_key] = scistats.skew(c_ells, bias=True)

            # * Read in EB data
            data_df = pd.read_csv(
                self.folder_path + self.masked_cl_out_dir + '/EB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                ell_key = int(label) - 2

                if mask_num == 4:
                    c_ells *= self.mask_f_sky

                # Compute the root-mean-square of the set of values.
                # As values can be negative, this is more accurate than the mean
                mean_val = np.sqrt(np.mean(c_ells ** 2))
                EB_rms[mask_key][ell_key] = mean_val

                # Compute variance
                EB_var[mask_key][ell_key] = np.var(c_ells ** 2) / (mean_val * mean_val)

                # Compute skew
                EB_skew[mask_key][ell_key] = scistats.skew(c_ells, bias=True)

        fig_size = [11, 6]

        # Create the colour-map from the f_sky values
        norm = mpl.colors.LogNorm(vmin=min(self.masks_f_sky), vmax=max(self.masks_f_sky))
        cmap = mpl.cm.ScalarMappable(norm=norm, cmap='plasma')
        cmap.set_array([])

        # Plot the EE mean
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 5):
            ax1.loglog(self.ells, EE_avg['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.loglog(self.ells, EE_avg['Mask5'], lw=2, c='cornflowerblue')

        # The expected EE signal for the shear given the convergence
        exp_EE = (self.ells + 2) * (self.ells - 1) / (self.ells * (self.ells + 1)) * self.c_ells['W4xW4'][2:]

        ax1.loglog(self.ells, exp_EE, lw=2, label=r'Expectation', c='cyan', ls='--')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\ell (\ell + 1) C_\ell^{\textrm{EE}} / 2 \pi$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig1.tight_layout()

        # Plot the EE ratio
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 5):
            ax1.semilogx(self.ells, EE_avg['Mask' + str(mask_num)] / exp_EE - 1, lw=2,
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.semilogx(self.ells, EE_avg['Mask5'] / exp_EE - 1, lw=2, c='cornflowerblue')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$C_\ell^{\textrm{EE}} / N(\ell) C_\ell^{\kappa \kappa} - 1$')
        ax1.set_title('Ratio of recovered EE signal to expected signal from convergence')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig1.tight_layout()

        # Plot the BB mean
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 4):
            ax1.loglog(self.ells, BB_avg['Mask' + str(mask_num)], lw=2,
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.loglog(self.ells, BB_avg['Mask5'], lw=2, c='cornflowerblue', label='Euclid')

        # Plot 1% of the expected signal
        ax1.loglog(self.ells, exp_EE * 0.01, lw=2, c='cyan', ls='--', label=r'$1\%$ of $EE$ signal')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\ell (\ell + 1) C_\ell^{\textrm{BB}} / 2 \pi$')
        ax1.set_title('$BB$ spectrum for masked maps only')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        ax1.legend()
        fig1.tight_layout()

        # Plot the EB RMS
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 4):
            ax1.loglog(self.ells, np.abs(EB_rms['Mask' + str(mask_num)]), lw=2,
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.loglog(self.ells, EB_rms['Mask5'], lw=2, c='cornflowerblue', label='Euclid')

        # Plot 1% of the expected signal
        ax1.loglog(self.ells, exp_EE * 0.01, lw=2, c='cyan', ls='--', label=r'$1\%$ of $EE$ signal')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\textrm{RMS}[\ell (\ell + 1) C_\ell^{\textrm{EB}} / 2 \pi]$')
        ax1.set_title('$EB$ spectrum for masked maps only')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        ax1.legend()
        fig1.tight_layout()

        # Plot the EE variance
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 5):
            ax1.loglog(self.ells, EE_var['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.loglog(self.ells, EE_var['Mask5'], lw=2, c='cornflowerblue')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\textrm{Var}[C_\ell^{\textrm{EE}} ] / \textrm{Avg}^2 [C_\ell^{\textrm{EE}}]$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig1.tight_layout()

        # Plot the BB variance
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 5):
            ax1.loglog(self.ells, BB_var['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.loglog(self.ells, BB_var['Mask5'], lw=2, c='cornflowerblue')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\textrm{Var}[C_\ell^{\textrm{BB}} ] / \textrm{Avg}^2 [C_\ell^{\textrm{BB}}]$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig1.tight_layout()

        # Plot EE, EB and BB signal for Euclid
        fig1, ax1 = plt.subplots(figsize=fig_size)

        ax1.loglog(self.ells, EE_avg['Mask5'], lw=2, c='cornflowerblue', label='$EE$')
        ax1.loglog(self.ells, 0.01 * EE_avg['Mask5'], lw=2, ls='--', c='cornflowerblue', label=r'$1\% \, EE$')
        ax1.loglog(self.ells, EB_rms['Mask5'], lw=2, c='orange', label='$EB$')
        ax1.loglog(self.ells, BB_avg['Mask5'], lw=2, c='hotpink', label='$BB$')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\ell (\ell + 1) C_\ell^{ij} / 2 \pi$')
        ax1.set_title('Summary of recovered power spectra for the Euclid mask')
        plt.legend()
        fig1.tight_layout()

        # Plot the ratio of EB to EE
        fig1, ax1 = plt.subplots(figsize=fig_size)

        for mask_num in range(1, 4):
            ax1.loglog(self.ells, EB_rms['Mask' + str(mask_num)] / EE_avg['Mask' + str(mask_num)], lw=2,
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.loglog(self.ells, EB_rms['Mask5'] / EE_avg['Mask5'], lw=2, c='cornflowerblue', label='Euclid')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$C_\ell^{EB} / C_\ell^{EE}$')
        ax1.set_title('Ratio of $EB$ to $EE$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        ax1.legend()
        fig1.tight_layout()

        # Plot the ratio of BB to EE
        fig1, ax1 = plt.subplots(figsize=fig_size)

        for mask_num in range(1, 4):
            ax1.loglog(self.ells, BB_avg['Mask' + str(mask_num)] / EE_avg['Mask' + str(mask_num)], lw=2,
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.loglog(self.ells, BB_avg['Mask5'] / EE_avg['Mask5'], lw=2, c='cornflowerblue', label='Euclid')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$C_\ell^{BB} / C_\ell^{EE}$')
        ax1.set_title('Ratio of $BB$ to $EE$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        ax1.legend()
        fig1.tight_layout()

        # Plot the ratio of BB and EB to EE
        fig1, ax1 = plt.subplots(figsize=fig_size)

        for mask_num in range(1, 4):
            ax1.loglog(self.ells, BB_avg['Mask' + str(mask_num)] / EE_avg['Mask' + str(mask_num)], lw=2,
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))
            ax1.loglog(self.ells, EB_rms['Mask' + str(mask_num)] / EE_avg['Mask' + str(mask_num)], lw=2, ls='--',
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.loglog(self.ells, BB_avg['Mask5'] / EE_avg['Mask5'], lw=2, c='cornflowerblue', label='$BB$')
        ax1.loglog(self.ells, EB_rms['Mask5'] / EE_avg['Mask5'], lw=2, ls='--', c='cornflowerblue', label='$EB$')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$C_\ell^{ij} / C_\ell^{EE}$')
        ax1.set_title('Ratio of $BB$ and $EB$ to $EE$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        ax1.legend()
        fig1.tight_layout()

        # Plot the ratio of BB to EB
        fig1, ax1 = plt.subplots(figsize=fig_size)

        for mask_num in range(1, 4):
            ax1.semilogx(self.ells, EB_rms['Mask' + str(mask_num)] / BB_avg['Mask' + str(mask_num)], lw=2,
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.semilogx(self.ells, EB_rms['Mask5'] / BB_avg['Mask5'] , lw=2, c='cornflowerblue', label='Euclid')

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\textrm{RMS}[C_\ell^{EB}] / \textrm{Avg}[C_\ell^{BB}]$')
        ax1.set_title('Ratio of $EB$ to $BB$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        ax1.legend()
        fig1.tight_layout()

        plt.show()

        # Plot the EB variance
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 5):
            ax1.loglog(self.ells, EB_var['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                       c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\textrm{Var}[C_\ell^{\textrm{EB}}] / \textrm{Avg}^2 [C_\ell^{\textrm{EB}}]$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig1.tight_layout()

        # Plot the EE skew
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 5):
            ax1.semilogx(self.ells, EE_skew['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\textrm{Skew}[C_\ell^{\textrm{EE}}]$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig1.tight_layout()

        # Plot the EB skew
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 5):
            ax1.semilogx(self.ells, EB_skew['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\textrm{Skew}[C_\ell^{\textrm{EB}}]$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig1.tight_layout()

        # Plot the BB skew
        fig1, ax1 = plt.subplots(figsize=fig_size)
        for mask_num in range(1, 5):
            ax1.semilogx(self.ells, BB_skew['Mask' + str(mask_num)], lw=2, label=r'Mask num ' + str(mask_num),
                         c=cmap.to_rgba(self.masks_f_sky[mask_num - 1]))

        ax1.set_xlabel(r'$\ell$')
        ax1.set_ylabel(r'$\textrm{Skew}[C_\ell^{\textrm{BB}}]$')
        fig1.colorbar(cmap, label=r'$f_\textrm{sky}$')
        fig1.tight_layout()

        plt.show()

    def compare_EB_lmax_Nside(self):
        num_maps = 4

        for mask_num in range(1, num_maps + 1):
            EE1 = []
            EB1 = []
            BB1 = []

            EE2 = []
            EB2 = []
            BB2 = []

            EE3 = []
            EB3 = []
            BB3 = []

            EE4 = []
            EB4 = []
            BB4 = []

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_1000_1024/Run1' + '/EE_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                EE1.append(np.mean(c_ells))

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_1000_1024/Run1' + '/EB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                EB1.append(np.sqrt(np.mean(c_ells ** 2)))

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_1000_1024/Run1' + '/BB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                BB1.append(np.mean(c_ells))

            """
            2048
            """

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_1000_2048/Run1' + '/EE_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                EE2.append(np.mean(c_ells))

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_1000_2048/Run1' + '/EB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                EB2.append(np.sqrt(np.mean(c_ells ** 2)))

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_1000_2048/Run1' + '/BB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                BB2.append(np.mean(c_ells))

            """
            2000, 1024
            """
            data_df = pd.read_csv('Data/Cuillin/EB_Shear_2000_1024/Run1' + '/EE_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                EE3.append(np.mean(c_ells))

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_2000_1024/Run1' + '/EB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                EB3.append(np.sqrt(np.mean(c_ells ** 2)))

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_2000_1024/Run1' + '/BB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                BB3.append(np.mean(c_ells))

            """
            2000, 2048
            """
            data_df = pd.read_csv('Data/Cuillin/EB_Shear_2000_2048/Run1' + '/EE_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                EE4.append(np.mean(c_ells))

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_2000_2048/Run1' + '/EB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                EB4.append(np.sqrt(np.mean(c_ells ** 2)))

            data_df = pd.read_csv('Data/Cuillin/EB_Shear_2000_2048/Run1' + '/BB_Cls_Mask' + str(mask_num) + '.csv')

            for label, c_ells in data_df.items():
                if mask_num is 4:
                    c_ells *= self.mask_f_sky

                BB4.append(np.mean(c_ells))

            ells1 = np.arange(2, 1000 + 1)
            ells2 = np.arange(2, 2000 + 1)

            # Plot the EE mean
            fig1, ax1 = plt.subplots(figsize=[11, 6])
            ax1.loglog(ells1, EE1, lw=2, label=r'$\ell = 1000, N = 1024$', c='tab:blue')
            ax1.loglog(ells1, EE2, lw=2, label=r'$\ell = 1000, N = 2048$', c='orange')
            ax1.loglog(ells2, EE3, lw=2, label=r'$\ell = 2000, N = 1024$', c='hotpink')
            ax1.loglog(ells2, EE4, lw=2, label=r'$\ell = 2000, N = 2048$', c='purple')

            # The expected EE signal for the shear given the convergence
            exp_EE = (self.ells + 2) * (self.ells - 1) / (self.ells * (self.ells + 1)) * self.c_ells['W4xW4'][2:]

            ax1.loglog(self.ells, exp_EE, lw=2, label=r'Expectation', c='cyan', ls='--')

            ax1.set_xlabel(r'$\ell$')
            ax1.set_ylabel(r'$\ell (\ell + 1) C_\ell^{\textrm{EE}} / 2 \pi$')
            plt.legend()
            fig1.tight_layout()

            fig1, ax1 = plt.subplots(figsize=[11, 6])
            ax1.loglog(ells1, EB1, lw=2, label=r'$\ell = 1000, N = 1024$', c='tab:blue')
            ax1.loglog(ells1, EB2, lw=2, label=r'$\ell = 1000, N = 2048$', c='orange')
            ax1.loglog(ells2, EB3, lw=2, label=r'$\ell = 2000, N = 1024$', c='hotpink')
            ax1.loglog(ells2, EB4, lw=2, label=r'$\ell = 2000, N = 2048$', c='purple')

            ax1.set_xlabel(r'$\ell$')
            ax1.set_ylabel(r'$\ell (\ell + 1) C_\ell^{\textrm{EB}} / 2 \pi$')
            plt.legend()
            fig1.tight_layout()

            fig1, ax1 = plt.subplots(figsize=[11, 6])
            ax1.loglog(ells1, BB1, lw=2, label=r'$\ell = 1000, N = 1024$', c='tab:blue')
            ax1.loglog(ells1, BB2, lw=2, label=r'$\ell = 1000, N = 2048$', c='orange')
            ax1.loglog(ells2, BB3, lw=2, label=r'$\ell = 2000, N = 1024$', c='hotpink')
            ax1.loglog(ells2, BB4, lw=2, label=r'$\ell = 2000, N = 2048$', c='purple')

            ax1.set_xlabel(r'$\ell$')
            ax1.set_ylabel(r'$\ell (\ell + 1) C_\ell^{\textrm{BB}} / 2 \pi$')
            plt.legend()
            fig1.tight_layout()
            plt.show()

    def pseudo_cl_recov(self):
        for mask_num in range(1, 4):
            mask = hp.read_map(self.folder_path + f'Masks/Mask{mask_num}.fits')
            f_sky = mask.sum() / mask.size

            mask_cls = hp.anafast(mask, lmax=2000)

            data_df = pd.read_csv(self.folder_path + self.masked_cl_out_dir + f'/TEB_Cls_{mask_num}.dat',
                                  sep=r'\s+')

            mask_matrix = camb.mathutils.pcl_coupling_matrix(mask_cls, lmax=1000)

            inv_mask_matrix = np.linalg.pinv(mask_matrix[2:, 2:])

            recov_TT = np.dot(inv_mask_matrix, data_df['TT'].to_numpy())

            fig1, ax1 = plt.subplots(figsize=(11, 6))
            ax1.loglog(self.ells, self.ells * (self.ells + 1) * data_df['TT'] / (2 * np.pi) / f_sky, lw=2,
                       c='tab:blue', label='No correction')
            ax1.loglog(self.ells, self.ells * (self.ells + 1) * recov_TT / (2 * np.pi), lw=2,
                       c='orange', label='Pseudo-Cl')
            ax1.loglog(self.ells, self.c_ells['W4xW4'][2:], lw=2, label=r'Expectation', c='cyan', ls='--')

            ax1.set_xlabel(r'$\ell$')
            ax1.set_ylabel(r'$\ell (\ell + 1) C_\ell^{\textrm{EE}} / 2 \pi$')
            fig1.legend()
            fig1.tight_layout()
            plt.show()

    def simple_likelihood_as_ns(self):
        """
        Function that implements a very simple likelihood code that tries to find the joint parameter constraints on
        A_s and n_s

        Returns:
            None
        """
        lmax = 2000
        ells = np.arange(2, lmax + 1)

        # The N_side parameter of the generated maps & masks
        n_side = 2048

        # Convert our N_side parameter to the number of pixels in the map
        n_pix = 12 * n_side ** 2

        intrinsic_gal_ellip = 0.21  # The standard deviation of the intrinsic galaxy ellipticity distribution

        avg_gal_den = 30  # This is the average surface galaxy density in [num gals / arc min^2]
        area_per_pix = 1.49E8 / n_pix  # This is the total area in arcmin^2 divided by the number of pixels
        num_gal_per_pix = avg_gal_den * area_per_pix

        # Generate random Gaussian noise that will be added to our maps
        random_noise = np.random.normal(loc=0, scale=intrinsic_gal_ellip / np.sqrt(num_gal_per_pix), size=n_pix)

        # Compute what the expected Cl VALUE (singular) is for the shape noise.
        theory_cl_noise = intrinsic_gal_ellip ** 2 / (avg_gal_den / (sciconst.arcminute ** 2))

        # Include the l(l+1) / 2pi factor in the Cl valeus
        theory_cl_noise = ells * (ells + 1) * theory_cl_noise / (2 * np.pi)

        # Optionally use no theory noise, and set all values to zero
        theory_cl_noise = np.zeros(lmax - 1)

        # Initiate a LCDM cosmology
        params = camb.CAMBparams()
        params.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112, mnu=0.06)
        params.InitPower.set_params(As=2.1E-9, ns=0.96)
        # params.set_for_lmax(lmax, lens_potential_accuracy=1)
        params.set_for_lmax(lmax, lens_potential_accuracy=0)

        params.Want_CMB = False
        params.WantCls = True
        params.NonLinear = camb.model.NonLinear_both

        # We want to evaluate the lensing power spectrum at z=2 only, for the moment.
        params.SourceWindows = [GaussianSourceWindow(redshift=2, source_type='lensing', sigma=0.05)]

        # The range of A_s values that we want to compute the power spectrum at
        a_s_vals = np.linspace(1.5e-9, 4e-9, 10)
        n_s_vals = np.linspace(0.7, 1.3, 10)

        # List which our Cl values will get stored into
        data_dict = {'As': [], 'ns': [], 'Cl': []}

        print('Evaluating CAMB on the grid now')

        # Go through each A_s value, compute the lensing power spec, and save to the list.
        for a_s_idx, a_s in enumerate(a_s_vals):
            for n_s_idx, n_s in enumerate(n_s_vals):
                params.InitPower.set_params(As=a_s, ns=n_s)

                results = camb.get_results(params)

                cls = results.get_source_cls_dict()

                data_dict['As'].append(a_s)
                data_dict['ns'].append(n_s)
                data_dict['Cl'].append(cls['W1xW1'][2:lmax + 1])

        print('Evaluate the splines now')

        # Now we want to spline the power spectrum at each ell to get this as a function of A_s.
        # This allows us to predict the power spectra at arbitrary ell.
        splines = []

        # The spline function that we will use to interpolate our 2D data
        spline_func = sciinterp.RectBivariateSpline

        for ell_idx in range(len(ells)):
            # Transform our list into a 1D array, which can then be reshaped into a 2D one
            c_l_arr = np.array([data_dict['Cl'][idx][ell_idx] for idx in range(len(a_s_vals) * len(n_s_vals))])

            splines.append(spline_func(a_s_vals, n_s_vals, np.reshape(c_l_arr, (len(a_s_vals), len(n_s_vals)))))

        # Read in the unmasked convergence map that was created with A_s = 2.25E-9 and m_nu = 0.12
        converg_map = hp.read_map(self.folder_path + 'KappaGammaMap-f2z2.fits', verbose=True, field=0)

        # Read in the Euclid mask
        euclid_mask = hp.read_map('./resources/Euclid_masks/Euclid-gal-mask-2048.fits', verbose=False).astype(
            np.bool)

        # Compute f_sky for Euclid mask
        f_sky = euclid_mask.sum() / euclid_mask.size

        # Create masked map from convergence map and Euclid mask
        converg_map_mask = hp.ma(converg_map)
        converg_map_mask.mask = np.logical_not(euclid_mask)

        # Add the random shape noise to our convergence map
        # converg_map += random_noise

        # Convert to C_ells for both unmasked and masked maps
        converg_cls = np.array(hp.anafast(converg_map, lmax=lmax)[2:])
        converg_cls = ells * (ells + 1) * converg_cls / (2 * np.pi)

        converg_cls_mask = np.array(hp.anafast(converg_map_mask, lmax=lmax)[2:])
        converg_cls_mask = ells * (ells + 1) * converg_cls_mask / (2 * np.pi)

        # Normalise masked Cl through 1/f_sky
        converg_cls_mask /= f_sky

        # Now we want to evaluate the likelihood for a range of A_s values to find the maximum-likelihood value
        a_s_vals = np.linspace(1.5e-9, 4e-9, 75)
        n_s_vals = np.linspace(0.7, 1.3, 75)

        # Likelihood for our unmasked, Euclid mask, and custom mask
        likelihoods = np.zeros(shape=(len(a_s_vals), len(n_s_vals)))
        likelihoods_mask = np.zeros(shape=(len(a_s_vals), len(n_s_vals)))

        print('Evaluating the likelihood now')

        for a_s_idx, a_s in enumerate(a_s_vals):
            for n_s_idx, n_s in enumerate(n_s_vals):
                # Compute theory values for this (a_s, n_s) combination
                cl_theory = np.array([splines[ell](a_s, n_s) for ell in range(len(ells))])
                cl_theory = cl_theory[:, 0, 0]

                # Compute the log-likelihood
                log_lik = -1 * np.sum((2 * ells + 1) * (np.log(cl_theory + theory_cl_noise) +
                                                        converg_cls / (cl_theory + theory_cl_noise)))
                likelihoods[a_s_idx, n_s_idx] = log_lik

                # Compute the log-likelihood for our masked map
                log_lik = -1 * np.sum((2 * ells + 1) * (np.log(cl_theory + theory_cl_noise) +
                                                        converg_cls_mask / (cl_theory + theory_cl_noise)))
                likelihoods_mask[a_s_idx, n_s_idx] = log_lik

        # Find what the maximum likelihood value is for both maps
        max_idx = np.unravel_index(np.argmax(likelihoods), likelihoods.shape)
        max_as = a_s_vals[max_idx[0]]
        max_ns = n_s_vals[max_idx[1]]

        max_idx_mask = np.unravel_index(np.argmax(likelihoods_mask), likelihoods_mask.shape)
        max_as_mask = a_s_vals[max_idx_mask[0]]
        max_ns_mask = n_s_vals[max_idx_mask[1]]

        print(f'Maximum likelihood point is at A_s = {max_as:.4e}, n_s = {max_ns:.4f} for the unmasked map')
        print(f'Maximum likelihood point is at A_s = {max_as_mask:.4e}, n_s = {max_ns_mask:.4f} for the masked map')

        # * Now save our 2D likelihoods for GetDist
        # To do so, we turn our 2D array into a Pandas DataFrame
        data_dict = {'weight': [], 'like': [], 'As': [], 'ns': []}

        for a_s_idx, a_s in enumerate(a_s_vals):
            for n_s_idx, n_s in enumerate(n_s_vals):
                data_dict['weight'].append(1.0)
                data_dict['like'].append(likelihoods[a_s_idx, n_s_idx])
                data_dict['As'].append(a_s)
                data_dict['ns'].append(n_s)

        # Turn the dictionary to a dataframe and save this to disk
        data_df = pd.DataFrame(data_dict)
        data_df.to_csv(self.folder_path + 'GetDist/values1.txt', sep=' ', index=False)

        # Remove data once done
        del data_dict
        del data_df

        import getdist
        from getdist import plots
        mpl.use('Qt5Agg')

        # Use GetDist to read in the samples that we've just saved - ignoring the header row
        samples = getdist.loadMCSamples(self.folder_path + 'GetDist/values1', settings={'ignore_rows': 1})

        # Create a triangle plot of the data
        g = plots.get_subplot_plotter()
        g.triangle_plot(samples, filled=True)

        # Save the plot
        g.export(self.folder_path + 'TrianglePlot.pdf')
        plt.show()

        likelihoods_as1 = []
        likelihoods_as2 = []
        likelihoods_as3 = []

        # As we've got fixed n_s here, just go through A_s
        for a_s_idx, a_s in enumerate(a_s_vals):
            # Compute theory values for n_s = 0.90
            cl_theory = np.array([splines[ell](a_s, 0.90) for ell in range(len(ells))])
            cl_theory = cl_theory[:, 0, 0]

            # Compute likelihoods for this n_s
            log_lik = -1 * np.sum((2 * ells + 1) * (np.log(cl_theory + theory_cl_noise) +
                                                    converg_cls / (cl_theory + theory_cl_noise)))
            likelihoods_as1.append(log_lik)

            # Compute theory values for n_s = 0.96
            cl_theory = np.array([splines[ell](a_s, 0.96) for ell in range(len(ells))])
            cl_theory = cl_theory[:, 0, 0]

            # Compute likelihoods for this n_s
            log_lik = -1 * np.sum((2 * ells + 1) * (np.log(cl_theory + theory_cl_noise) +
                                                    converg_cls / (cl_theory + theory_cl_noise)))
            likelihoods_as2.append(log_lik)

            # Compute theory values for n_s = 1.05
            cl_theory = np.array([splines[ell](a_s, 1.05) for ell in range(len(ells))])
            cl_theory = cl_theory[:, 0, 0]

            # Compute likelihoods for this n_s
            log_lik = -1 * np.sum((2 * ells + 1) * (np.log(cl_theory + theory_cl_noise) +
                                                    converg_cls / (cl_theory + theory_cl_noise)))
            likelihoods_as3.append(log_lik)

        # Compute percentiles
        percent_68 = np.percentile(likelihoods.reshape(likelihoods.size), 32)
        percent_95 = np.percentile(likelihoods.reshape(likelihoods.size), 5)

        print(f'68-th percentile: {percent_68:.4e}')
        print(f'95-th percentile: {percent_95:.4e}')

        # Create a 1D histogram of the likelihood values, and draw on the percentiles too
        plt.figure()
        sns.histplot(likelihoods.reshape(likelihoods.size), label='Unmasked data', color='tab:blue')
        sns.histplot(likelihoods_mask.reshape(likelihoods_mask.size), label='Masked data', color='tab:pink')

        plt.axvline(percent_68, c='cornflowerblue', ls='--', lw=2, label='68-th percentile')
        plt.axvline(percent_95, c='orange', ls='--', lw=2, label='95-th percentile')
        plt.xlabel(r'$\ln \mathcal{L}$')
        plt.legend()
        plt.tight_layout()

        # Plot the contour of the likelihood values, along with the percentiles
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(13, 10), sharey=True)

        min_val = np.min([likelihoods, likelihoods_mask])
        max_val = np.max([likelihoods, likelihoods_mask])

        # Ax1 is for the unmasked values
        cont1 = ax1.contourf(a_s_vals, n_s_vals, likelihoods, levels=250, cmap='viridis',
                             vmin=min_val, vmax=max_val)

        c1 = ax1.contour(a_s_vals, n_s_vals, likelihoods, levels=[percent_68], colors='cornflowerblue', linestyles='--',
                         linewidth=2)
        c2 = ax1.contour(a_s_vals, n_s_vals, likelihoods, levels=[percent_95], colors='orange', linestyles='--',
                         linewidth=2)

        c1.collections[0].set_label('68-th percentile')
        c2.collections[0].set_label('95-th percentile')

        ax1.plot(max_as, max_ns, 'x', c='hotpink', lw=4, markersize=10, label='Maximum likelihood')
        ax1.plot(2.1E-9, 0.96, 'x', c='purple', lw=4, markersize=10, label='True value')

        # Ax2 is for the masked values
        cont2 = ax2.contourf(a_s_vals, n_s_vals, likelihoods_mask, levels=250, cmap='viridis',
                             vmin=min_val, vmax=max_val)

        ax2.plot(max_as_mask, max_ns_mask, 'x', c='hotpink', lw=4, markersize=10)
        ax2.plot(2.1E-9, 0.96, 'x', c='purple', lw=4, markersize=10)

        fig.subplots_adjust(right=0.85)
        cbar_ax = fig.add_axes([0.9, 0.15, 0.05, 0.7])
        fig.colorbar(cont1, label=r'$\ln \mathcal{L}$', cax=cbar_ax)

        ax1.set_xlabel(r'$A_\textrm{s}$')
        ax2.set_xlabel(r'$A_\textrm{s}$')
        ax1.set_ylabel(r'$n_\textrm{s}$')

        ax1.set_title('Unmasked map')
        ax2.set_title('Masked map')

        fig.suptitle(r'Log-likelihoods in the $A_\textrm{s} - n_\textrm{s}$ plane')
        fig.legend()

        # * Plot the difference in likelihood between masked and unmasked maps
        fig, ax = plt.subplots(figsize=(11, 6))
        cont = ax.contourf(a_s_vals, n_s_vals, likelihoods - likelihoods_mask,
                           levels=250, cmap=make_planck_colour_map())
        ax.set_xlabel(r'$A_s$')
        ax.set_ylabel(r'$n_s$')
        ax.set_title(r'Difference in log-likelihood between unmasked and masked maps')
        fig.colorbar(cont, label=r'$\ln \mathcal{L} - \ln \mathcal{L}_\textrm{mask}$', aspect=15)

        # Now plot the 1D-likelihoods at fixed n_s
        fig, ax = plt.subplots(figsize=(11, 6))
        ax.scatter(a_s_vals, likelihoods_as1, c='orange', label=r'$n_\textrm{s} = 0.90$')
        ax.scatter(a_s_vals, likelihoods_as2, c='tab:blue', label=r'$n_\textrm{s} = 0.96$')
        ax.scatter(a_s_vals, likelihoods_as3, c='hotpink', label=r'$n_\textrm{s} = 1.05$')
        ax.set_xlabel(r'$A_s$')
        ax.set_ylabel(r'$\ln \mathcal{L}$')
        ax.set_title(r'Log-likelihood as a function of $A_\textrm{s}$ for fixed $n_\textrm{s}$')
        plt.legend()
        plt.tight_layout()

        plt.show()
