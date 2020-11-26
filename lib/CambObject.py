import os
import subprocess
import random
import pathlib
import time
import ctypes
import numpy as np
from scipy import stats as scistats
import pandas as pd
import camb
import healpy as hp
import matplotlib.pyplot as plt
import seaborn as sns

from .flask_scripts.prepCambInput import split_files
from .flask_scripts.camb2info import XavierShift
from .Planck_colourmap import make_planck_colour_map


sns.set(font_scale=1.75, rc={'text.usetex': True})


class CambObject:

    def __init__(self, folder_path, lmax, non_linear=True):
        """
        General class that contains functions necessary to first compute the lensing power spectra for a given
        cosmology, save this power spectrum to the disk, input this spectra into Flask, run Flask multiple times
        to obtain multiple realisations of the Cl coefficients, and plotting the results.

        Args:
            folder_path (str): String that identifies where output data should be stored under the ./Data/ sub-directory
            lmax (int): Integer specifying the maximum l value that the lensing power spectrum should be computed to.
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
        return mean_dens * (z / z_0) ** 2 * np.exp(-(z / z_0) ** (3 / 2))

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

    def write_flask_config_file(self, n_side=2048):
        """
        Function that writes a Flask configuration file to the disk

        Args:
            n_side (int): The N_side parameter that is used to construct the maps in Flask

        Returns:
            None

        """

        # Do we want Flask to run with additional inputs that use the expected galaxy number densities to add noise
        # the the lensing signal?
        with_galaxy_counts = self.galaxy_dens

        filename = 'FlaskInput.ini'
        filename = self.folder_path + filename

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
        file.write('STARMASK: \t 0 \n\n')  # TODO: change this when have a mask!

        file.write('\n## Multipole information ##\n\n')
        file.write('EXTRAP_DIPOLE: \t 0 \n')
        file.write('LRANGE: \t 2 ' + str(self.ell_max) + '\n')
        file.write('CROP_CL: \t 0 \n')
        file.write('SHEAR_LMAX: \t' + str(self.ell_max) + '\n')
        file.write('NSIDE: \t ' + str(n_side) + ' \n')
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
            file.write('EXIT_AT: \t ELLIPFITS_PREFIX \n')
        else:
            file.write('EXIT_AT: \t RECOVCLS_OUT \n')

        file.write('FITS2TGA: \t 0 \n')
        file.write('USE_UNSEEN: \t 1 \n')
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
        file.write('RECOVCLS_OUT: \t Output-Cl.dat \n')

        # Here, we have the shear values output the convergence and shear values at each redshift bin of the lensing
        # field, but with noise only due to cosmic variance, NOT due to the mask or Poisson galaxy sampling.
        file.write('SHEAR_ALM_PREFIX: \t 0 \n')
        file.write('SHEAR_FITS_PREFIX: \t 0 \n')  # output is (kappa, gamma1, gamma2)
        file.write('SHEAR_MAP_OUT: \t 0 \n')

        # This is the map data of each field AFTER applying cosmic variance, mask, and Poisson galaxy sampling.
        # i.e. this is the most "noisiest" output, but is the closest to the observed values and what we need to model
        if with_galaxy_counts:
            file.write('MAPWERFITS_PREFIX: \t Output-poisson-map- \n')
        else:
            file.write('MAPWERFITS_PREFIX: \t 0 \n')
        file.write('MAPWER_OUT: \t 0 \n')

        # Here outputs galaxy ellipticities values at each coordinate that includes galaxy noise
        file.write('ELLIPFITS_PREFIX: \t Output-Ellip-map- \n')
        file.write('ELLIP_MAP_OUT: \t 0 \n')

        file.write('CATALOG_OUT: \t 0 \n')
        file.write('\nCATALOG_COLS: \t theta phi z kappa gamma1 gamma2 ellip1 ellip2\n')

        file.close()

    def run_flask(self, use_rnd=None):
        """
        Function that runs the Flask executable on an already existing Flask input .ini file

        Args:
            use_rnd (bool): Flag which, when enabled, generates a new random number to be used as Flask's random
                            number seed - which allows the same input .ini file to generate multiple different spectra

        Returns:
            None, all data is written to the disk by Flask
        """
        # Check that the Flask executable path exists before running it
        if self.flask_executable is None:
            raise RuntimeError('Please first set the location of the Flask executable before running it!')

        if not os.path.isfile(self.folder_path + 'FlaskInput.ini'):
            print('The Flask input file has not been previously generated, making one now')
            self.write_flask_config_file()

        print('\nRunning Flask')

        # Time how long it takes Flask to run
        start_time = time.time()

        # Execute Flask as a subprocess run from the shell.
        command = subprocess.run(str(self.flask_executable) + ' FlaskInput.ini ' +
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
            output_file.close()

        # Else, there was an error somewhere, print this error and exit the program
        else:
            print('Flask did not run successfully, hopefully the error is contained below')
            print(command.stdout)
            print(command.stderr)

            raise RuntimeError('Flask did not run successfully! :(')

    def multiple_run_flask(self, num_runs):
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

        # Run Flask the given number of times to generate a power spectra each time
        for run_num in range(num_runs):
            self.run_flask(use_rnd=True)

            # TODO: change these labels to deal with galaxy counts too!
            # Read in the new power spectrum
            cl_df = pd.read_csv(self.folder_path + 'Output-Cl.dat', header=None,
                                names=['ell', 'Cl-f1z1f1z1', 'Cl-f1z1f1z2', 'Cl-f1z2f1z2'], sep=r'\s+', skiprows=1)

            # Strip any spaces out of the column names
            cl_df.columns = cl_df.columns.str.lstrip()

            # Get the ell values as a numpy array
            ells = cl_df['ell'].to_numpy()

            # Create a new data-frame that stores the raw Cl values (i.e. without ell(ell+1) normalisation)
            raw_cl_df = pd.DataFrame({'Cl': cl_df['Cl-f1z2f1z2'].to_numpy()}, index=ells)

            # Get the c_ell values normalised by ell(ell+1)/2pi, as usual
            c_ells = cl_df['Cl-f1z2f1z2'].to_numpy() * self.ells * (self.ells + 1) / (2 * np.pi)

            # Create a new dummy dataframe with out C_ell values, with an index determined by the ell values
            cl_df = pd.DataFrame({'Cl': c_ells}, index=ells)

            # Append the transpose of the current data-frames to the lists
            raw_data_dfs.append(raw_cl_df.transpose())
            data_dfs.append(cl_df.transpose())

        # Concatenate the all the data frames together into a single output data frame, for both data-sets
        raw_data_df = pd.concat(raw_data_dfs, ignore_index=True)
        data_df = pd.concat(data_dfs, ignore_index=True)

        # Prints timing statistics
        print('The total time for the runs was {num:.3f} seconds, with an average of {numpersec:.3f} seconds/run'
              .format(num=time.time() - start_time, numpersec=(time.time() - start_time) / num_runs))

        # Save the data which can then be read in later
        raw_data_df.to_csv(self.folder_path + 'AggregateRawCls.csv', index=False)
        data_df.to_csv(self.folder_path + 'AggregateCls.csv', index=False)

        # Hand off to the plotting function which plots the output data
        self.plot_multiple_run_data()

    def plot_multiple_run_data(self):
        """
        Function that plots data that has been calculated for multiple runs of Flask

        Returns:
            None
        """
        # Import the previously saved data
        data_df = pd.read_csv(self.folder_path + 'AggregateRawCls.csv')
        data_df_kde = pd.read_csv(self.folder_path + 'AggregateCls4.csv')

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

        # Get the keys of the data frame
        keys = data_df.keys()

        # Dictionary which the Pearson data will be stored into
        pearson_data = {'x': [], 'y': [], 'r': []}

        covariance_data = {'x': [], 'y': [], 'raw_C': [], 'norm_C': []}

        # Go through each Cl combination and work out the Pearson correlation coefficient
        for ell1, c_ells1 in data_df.items():
            # Only go up to ell1 of 125, otherwise figure too crowded
            if int(ell1) > 126:
                continue

            # Only want to compute samples where ell2 < ell1
            for ell2, c_ells2 in data_df.items():
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
                # Go through each ell1 and ell2 list and compute the sum at that point
                for val1, val2 in zip(c_ells1.to_numpy(), c_ells2.to_numpy()):
                    val += (val1 - self.raw_c_ells['W2xW2'][int(ell1)]) * \
                           (val2 - self.raw_c_ells['W2xW2'][int(ell2)])

                # Normalise the sum to the number of data points
                val /= len(c_ells1)

                # Store the data in the covariance data dictionary
                covariance_data['x'].append(int(ell1))
                covariance_data['y'].append(int(ell2))
                covariance_data['raw_C'].append(val)

                # Now normalise through the intrinsic standard deviation at each ell value
                val /= self.raw_c_ells['W2xW2'][int(ell1)] * np.sqrt(2 / (2 * int(ell1) + 1))
                val /= self.raw_c_ells['W2xW2'][int(ell2)] * np.sqrt(2 / (2 * int(ell2) + 1))

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
            ax = sns.heatmap(data=norm_covariance_df, square=True, cmap="jet",
                             cbar_kws={'label': r'Covariance matrix $C_{ij} / \sigma_i \sigma_j$'})
            ax.set_xlabel(r'$\ell_i$')
            ax.set_ylabel(r'$\ell_j$')
            ax.set_title('Normalised covariance matrix')
            plt.tight_layout()

            plt.show()

        # In order to plot a KDE, we need to increase the size of the values in order to make them order-unity
        for idx in range(len(keys)):
            data_df_kde[data_df_kde.keys()[idx]] = (data_df_kde[data_df_kde.keys()[idx]] * 1E7)

        grid2 = sns.PairGrid(data_df_kde, vars=[keys[0], keys[2], keys[8], keys[23], keys[98], keys[498], keys[1498]])
        grid2.map_diag(sns.histplot)
        grid2.map_lower(sns.kdeplot, shade=True, levels=4)
        grid2.map_upper(sns.histplot)
        grid2.tight_layout()
        plt.show(block=False)

        # * Print the raw data for each Cl dataset - very messy but gets the point across
        plt.figure(figsize=(13, 7))
        for cl in data_df.itertuples(index=False):
            plt.loglog(self.ells, cl, alpha=0.5, linewidth=0.75)

        plt.loglog(self.ells, mean_cls, lw=2, color='blue', label=r'Average $C_\ell$')

        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')

        plt.legend()
        plt.tight_layout()

        # * Plot the deviation of the mean_cls vs input cls from CAMB
        plt.figure(figsize=(13, 7))
        plt.semilogx(self.ells, (np.array(mean_cls) / self.raw_c_ells['W2xW2'][2:]) - 1, lw=2, color='purple')

        plt.title(r'Deviation of the average $C_\ell$ with respect to the input $C_\ell$ values')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$ (C_\ell^\textrm{Avg} - C_\ell^\textrm{In} ) / C_\ell^\textrm{In} $')

        plt.tight_layout()

        # Also plot the skew
        plt.figure(figsize=(13, 7))

        # Plot the skew
        plt.semilogx(self.ells, skew_cls, lw=1, color='b', label='Fully numerical')
        plt.semilogx(self.ells, skew_cls2, lw=1, color='k', label='With theory variance')

        # Plot the expected skew for a Gamma distribution
        k = (2 * self.ells + 1) / 2
        plt.semilogx(self.ells, 2 / np.sqrt(k), color='purple', lw=2, label=r'$\Gamma$ function prediction')

        plt.title('Skew')
        plt.xlabel(r'$\ell$')
        plt.legend()
        plt.tight_layout()

        # And kurtosis
        plt.figure(figsize=(13, 7))
        plt.semilogx(self.ells, kurt_cls, lw=1, color='b', label='Fully numerical')
        plt.semilogx(self.ells, kurt_cls2, lw=1, color='k', label='With theory variance')
        plt.semilogx(self.ells, 6 / k, color='purple', lw=2, label=r'$\Gamma$ function prediction')
        plt.title('Kurtosis')
        plt.xlabel(r'$\ell$')
        plt.legend()
        plt.tight_layout()

        # Now plot the variance of the Cl's with the expected values from the Cosmic variance
        plt.figure(figsize=(13, 7))
        plt.loglog(self.ells, var_cls, lw=2, color='purple', label='Data')
        plt.loglog(self.ells, 2 / (2 * self.ells + 1), color='blue', label='Cosmic variance')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\textrm{Var}[C_\ell] / \textrm{Avg}[C_\ell]^2$')
        plt.legend()
        plt.tight_layout()
        plt.show()

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
        data_df = pd.read_csv(self.folder_path + 'AggregateCls1.csv')

        # Normalise each value by 1E7 to get working KDE plots
        for idx in range(len(data_df.keys())):
            data_df[data_df.keys()[idx]] = (data_df[data_df.keys()[idx]] * 1E7)

        data_frames = []

        # Go through each ell and extract the Cl data
        for label, c_ells in data_df.items():
            # Subject to the condition where ell <= 25
            if int(label) > 25:
                continue
            tmp_df = pd.DataFrame({'ell': label, 'Cl': (c_ells - np.mean(c_ells)) / np.mean(c_ells)})
            data_frames.append(tmp_df)

        data = pd.concat(data_frames, ignore_index=True)

        # Resets the Seaborn pallet to have a white figure background
        sns.set_theme(style="white", rc={"axes.facecolor": (0, 0, 0, 0)})

        # Create the colour map for our ell values
        pal = sns.cubehelix_palette(len(data_frames), rot=-0.2, light=0.75, dark=0.05, gamma=1.25)

        # Initialize the FacetGrid object
        g = sns.FacetGrid(data, row="ell", hue="ell", aspect=20, height=0.5, palette=pal)

        # Draw the densities in a few steps
        g.map(sns.kdeplot, "Cl", bw_adjust=0.5, clip_on=False, fill=True, alpha=1, linewidth=1.5)
        g.map(sns.kdeplot, "Cl", clip_on=False, color="w", lw=2, bw_adjust=0.5)
        g.map(plt.axhline, y=0, lw=2, clip_on=False)
        g.map(plt.axvline, x=0, linestyle='--', color='w', lw=1, alpha=0.8, clip_on=False)

        # Define and use a simple function to label the plot in axes coordinates
        def label_func(x, color, label):
            ax = plt.gca()
            ax.text(0, .2, label, fontweight="bold", color=color,
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
        g.map(plot_gamma, "Cl")

        # Set the subplots to overlap
        g.fig.subplots_adjust(hspace=-0.45)

        # Add title to figure and adjust its position
        g.fig.subplots_adjust(top=0.995)
        g.fig.suptitle(r'Distribution of the $C_\ell$ values', fontsize=18)
        plt.xlabel(r'$[C_\ell - \bar{C}_\ell] / \bar{C}_\ell$', fontsize=16)
        plt.xlim(left=-1, right=1)

        # Remove axes details that don't play well with overlap
        g.set_titles("")
        g.set(yticks=[])
        g.despine(bottom=True, left=True)
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

    def use_cpp_map_to_alm(self):
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

    def experimenting_with_masks(self):
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
        plt.show(block=False)

        n_side = 2048

        n_pix = 12 * n_side ** 2
        random_map1 = np.random.normal(loc=0, scale=1, size=n_pix)
        random_map2 = np.random.normal(loc=0, scale=np.sqrt(10), size=n_pix)
        hp.mollview(random_map1, title="Gasussian random noise in each pixel", norm="hist", cmap='viridis')
        hp.graticule(verbose=False)

        cl1 = hp.anafast(random_map1, lmax=2000)
        cl2 = hp.anafast(random_map2, lmax=2000)
        ell = np.arange(2, 2000 + 1)

        plt.figure(figsize=(13, 7))
        plt.loglog(ell, ell * (ell + 1) * cl1[2:] / (2 * np.pi), label=r'$\sigma^2 = 1$')
        plt.loglog(ell, ell * (ell + 1) * cl2[2:] / (2 * np.pi), label=r'$\sigma^2 = 10$')
        plt.xlabel(r"$\ell$")
        plt.ylabel(r"$\ell(\ell+1)C_{\ell} / 2 \pi$")
        plt.title('Power spectrum of Gaussian random noise')

        plt.legend()
        plt.tight_layout()
        plt.show()

        intrinsic_gal_ellip = 0.21  # The standard deviation of the intrinsic galaxy ellipticity distirbution
        avg_gal_den = 30  # This is the average surface galaxy density in [num gals / arc min^2]
        area_per_pix = 1.49E8 / n_pix  # This is the total area in arcmin^2 divided by the number of pixels
        num_gal_per_pix = avg_gal_den * area_per_pix

        sigma_noise = num_gal_per_pix * intrinsic_gal_ellip

        print('sigma_noise is: ', sigma_noise)

        new_map = np.zeros(shape=n_pix)

        for i in range(len(new_map)):
            new_map[i] = np.random.normal(loc=0, scale=sigma_noise)

        print(new_map)

        hp.mollview(new_map)
        hp.graticule(verbose=False)
        plt.show()

        noise_cl = hp.anafast(random_map1, lmax=2000)
        plt.loglog(ell, ell * (ell + 1) * noise_cl[2:] / (2 * np.pi), label=r'$\sigma^2 = 1$')
        plt.show()

        import sys
        sys.exit()
