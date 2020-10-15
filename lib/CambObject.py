import os
import subprocess
import pathlib
import time
import ctypes
import numpy as np
import pandas as pd
import camb
import healpy as hp
import matplotlib.pyplot as plt

from .flask_scripts.prepCambInput import split_files
from .flask_scripts.camb2info import XavierShift
from .Planck_colourmap import make_planck_colour_map


class CambObject:

    def __init__(self, folder_path, lmax, non_linear=True):
        # TODO: a rework of data sorting in folders
        self.folder_path = 'Data/' + folder_path + '/'
        if not os.path.isdir(self.folder_path):
            os.makedirs(self.folder_path)

        self.ell_max = lmax
        self.non_linear = non_linear

        self.ells = np.arange(2, lmax + 1)

        self.params = camb.CAMBparams()
        self.params.set_cosmology(H0=70, ombh2=0.0226, omch2=0.112)
        self.params.InitPower.set_params(As=2.1E-9, ns=0.96)
        self.params.set_for_lmax(self.ell_max, lens_potential_accuracy=1)

        self.params.Want_CMB = True

        if self.non_linear:
            self.params.NonLinear = camb.model.NonLinear_both
            self.params.set_nonlinear_lensing(True)

        else:
            self.params.set_nonlinear_lensing(False)

        self.window_functions = None
        self.num_redshift_bins = None
        self.results = None
        self.c_ells = None
        self.raw_c_ells = None
        self.camb_output_filename = None

        self.flask_executable = None
        self.flask_executed = False

        # A list that holds a list for each redshift bin of the Cl's computed by me
        self.my_cls = []
        self.my_cls_cpp = []

    def set_flask_executable(self, path):
        """
        Function to set the executable path for the FLASK software

        :param path: The path of the parameter

        :return: None
        """

        # Convert the given path to an absolute path
        path = pathlib.Path(path).expanduser().resolve()

        # First check that the file exists
        if not path.is_file():
            raise RuntimeError('The given path for FLASK is not correct! Please check.')

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
        if not isinstance(window_functions, list):
            raise RuntimeError('The argument "window_functions" needs to be a list of window functions')

        self.window_functions = window_functions
        self.num_redshift_bins = len(window_functions)

        self.params.SourceWindows = [window_func.construct_camb_instance() for window_func in window_functions]

    def compute_c_ells(self):
        start_time = time.time()
        print('Running CAMB')

        self.results = camb.get_results(self.params)
        self.c_ells = self.results.get_source_cls_dict(lmax=self.ell_max)
        self.raw_c_ells = self.results.get_cmb_unlensed_scalar_array_dict(lmax=self.ell_max, raw_cl=True,
                                                                          CMB_unit='muK')

        print('CAMB finished in {num:.2f} seconds'.format(num=time.time() - start_time))

    def get_c_ells_dict(self, key=None):
        if self.c_ells is None:
            self.compute_c_ells()

        if key is None:
            return self.c_ells

        else:
            return self.c_ells[key]

    def output_c_ells(self):
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
        split_files(self.camb_output_filename,
                    self.camb_output_filename.split('.')[0] + '-',
                    self.num_redshift_bins)

    def plot_1x2_map_power_spectrum(self, key, nside=512):
        if key not in self.raw_c_ells.keys():
            raise RuntimeError('The selected key "' + str(key) + '" for plotting is not in the computed array of: '
                               + str(self.raw_c_ells.keys()))

        # Obtain the Cl values using the provided key for the results dictionary. Only plot ell >= 2
        cl = self.raw_c_ells[key][2:]

        # Obtain the Planck colour map, which we will use in plots
        planck_cmap = make_planck_colour_map()

        # Use Healpy to plot the Cl's as a map
        map1 = hp.sphtfunc.synfast(cl, nside=nside, new=True, verbose=False)
        hp.visufunc.mollview(map1, cmap=planck_cmap, title='Map of power spectrum for ' + str(key))
        hp.graticule(verbose=False, alpha=0.6)
        plt.show()

        # Now convert our map to Cl's for comparison with the original Cl's
        cl_from_map = hp.sphtfunc.anafast(map1, lmax=self.ell_max)

        # Now, we plot both sets of Cl's
        plt.figure(figsize=(13, 7))
        plt.loglog(np.arange(2, self.ell_max + 1), cl, label=r'$C_\ell$ input', lw=2.5)
        plt.loglog(np.arange(2, self.ell_max + 1), cl_from_map[2:self.ell_max + 1], label=r'$C_\ell$ from map')
        plt.title(r'Comparison of $C_{\ell}$')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$C_{\ell}$')
        plt.legend()
        plt.tight_layout()
        plt.show()

        ells = np.arange(2, self.ell_max + 1)

        # Plot both of the Cl's, but now with the scaling of ell(ell + 1)/2pi in place
        plt.figure(figsize=(13, 7))
        plt.loglog(np.arange(2, self.ell_max + 1), cl * ells * (ells + 1) / (2 * np.pi), label=r'$C_\ell$ input',
                   lw=2.5)
        plt.loglog(np.arange(2, self.ell_max + 1), cl_from_map[2:self.ell_max + 1] * ells * (ells + 1) / (2 * np.pi),
                   label=r'$C_\ell$ from map')
        plt.title(r'Comparison of $C_{\ell}$')
        plt.xlabel(r'$\ell$')
        plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
        plt.legend()
        plt.tight_layout()
        plt.show()

    def create_fields_info(self):
        filename = 'fields_info.ini'

        file = open(self.folder_path + '/' + filename, 'w')

        file.write('# File auto-generated which contains information about the redshift distributions\n\n')
        file.write('# Field number, z bin number, mean, shift, field type (1:galaxy, 2:shear), zmin, zmax\n\n')

        for index, window in enumerate(self.window_functions, start=1):
            file.write('\t1 \t ' + str(index) + '\t0.0000\t' + '{:.4f}'.format(XavierShift(window.redshift))
                       + '\t2\t' + '{:.4f}'.format(window.redshift - window.sigma)
                       + '\t' + '{:.4f}'.format(window.redshift + window.sigma) + '\n')

        file.close()

    def write_flask_config_file(self):
        filename = 'FlaskInput.ini'
        filename = self.folder_path + filename

        file = open(filename, 'w')

        file.write('# This is a config file auto-generated by my code\n')
        file.write('# This file contains all of the cosmological, input and output data for the galaxy lensing\n\n\n')

        file.write('## Simulation basics ##\n\n')
        file.write('DIST: \t GAUSSIAN \n')
        file.write('RNDSEED: \t 315 \n')  # TODO: change this randomly?!
        file.write('POISSON: \t 1 \n')

        file.write('\n## Cosmology ##\n\n')
        file.write('OMEGA_m: \t 0.3 \n')
        file.write('OMEGA_L: \t 0.7 \n')
        file.write('W_de: \t -1.0 \n\n')
        file.write('ELLIP_SIGMA: \t 0.11 \n')
        file.write('GALDENSITY: \t 30 \n\n')

        file.write('\n## Input data ## \n\n')
        file.write('FIELDS_INFO: \t fields_info.ini \n')
        file.write('CHOL_IN_PREFIX: \t 0 \n')
        file.write('CL_PREFIX: \t PowerSpecCambOut- \n')
        file.write('ALLOW_MISS_CL: \t 0 \n')
        file.write('SCALE_CLS: \t 1.0 \n')
        file.write('WINFUNC_SIGMA: \t -1 \n')
        file.write('APPLY_PIXWIN: \t 0 \n')  # TODO: changed this
        file.write('SUPPRESS_L: \t -1 \n')
        file.write('SUP_INDEX: \t -1 \n\n')

        file.write('\n## Survey selection functions ##\n\n')
        file.write('SELEC_SEPARABLE: \t 1 \n')
        file.write('SELEC_PREFIX: \t 0 \n')
        file.write('SELEC_Z_PREFIX: \t example-z-selection- \n')
        file.write('SELEC_SCALE: \t 1 \n')
        file.write('SELEC_TYPE: \t 0 \n')
        file.write('STARMASK: \t 0 \n\n')

        file.write('\n## Multipole information ##\n\n')
        file.write('EXTRAP_DIPOLE: \t 0 \n')
        file.write('LRANGE: \t 2 ' + str(self.ell_max) + '\n')
        file.write('CROP_CL: \t 0 \n')
        file.write('SHEAR_LMAX: \t' + str(self.ell_max) + '\n')
        file.write('NSIDE: \t 2048 \n')  # TODO: change this?!
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
        file.write('EXIT_AT: \t 0 \n')
        file.write('FITS2TGA: \t 0 \n')
        file.write('USE_UNSEEN: \t 1 \n')
        file.write('LRANGE_OUT: \t 2 ' + str(self.ell_max) + '\n')
        file.write('MMAX_OUT: \t -1 \n')
        file.write('ANGULAR_COORD: \t 0 \n')
        file.write('DENS2KAPPA: \t 0 \n')

        file.write('FLIST_OUT: \t 0 \n')
        file.write('SMOOTH_CL_PREFIX: \t Output-smoothed-cl- \n')
        file.write('XIOUT_PREFIX: \t 0 \n')
        file.write('GXIOUT_PREFIX: \t 0 \n')
        file.write('GCLOUT_PREFIX: \t 0 \n')
        file.write('COVL_PREFIX: \t 0 \n')
        file.write('REG_COVL_PREFIX: \t 0 \n')
        file.write('REG_CL_PREFIX: \t Output-Reg-Cl- \n')  # Changed this
        file.write('CHOLESKY_PREFIX: \t 0 \n')
        file.write('AUXALM_OUT: \t 0 \n')  # TODO: change this?!
        file.write('RECOVAUXCLS_OUT: \t 0 \n')
        file.write('AUXMAP_OUT: \t 0 \n')
        file.write('DENS2KAPPA_STAT: \t 0 \n')
        file.write('MAP_OUT: \t 0 \n')
        file.write('RECOVALM_OUT: \t 0 \n')
        file.write('RECOVCLS_OUT: \t 0 \n')
        file.write('SHEAR_ALM_PREFIX: \t Output-shear-alm- \n')

        file.write('MAPFITS_PREFIX: \t Output-converg-map- \n')
        file.write('SHEAR_FITS_PREFIX: \t Output-shear-map- \n')
        file.write('MAPWERFITS_PREFIX: \t Output-poisson-map- \n')
        file.write('ELLIPFITS_PREFIX: \t Output-ellip-map- \n')

        file.write('SHEAR_MAP_OUT: \t 0 \n')
        file.write('MAPWER_OUT: \t 0 \n')
        file.write('ELLIP_MAP_OUT: \t 0 \n')
        file.write('CATALOG_OUT: \t 0 \n')

        file.write('\nCATALOG_COLS: \t theta phi z kappa gamma1 gamma2 ellip1 ellip2\n')

        file.close()

    def run_flask(self):
        if self.flask_executable is None:
            raise RuntimeError('Please first set the location of the Flask executable before running it!')

        print('Running Flask')

        start_time = time.time()

        # TODO: consider using a random, random number here for the RNDSEED value as this may introduce differences
        #  between the runs?

        command = subprocess.run([str(self.flask_executable), 'FlaskInput.ini'], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, universal_newlines=True, cwd=self.folder_path)

        if command.returncode == 0:
            print('Flask ran successfully :)')
            print('Flask took {num:.2f} seconds to run'.format(num=time.time() - start_time))
            self.flask_executed = True

            # Write Flask's output to a file
            output_file = open(self.folder_path + 'FlaskOutput.txt', 'w')
            output_file.write(command.stdout)
            output_file.close()

        else:
            print('Flask did not run successfully, hopefully the error is contained below')
            print(command.stdout)
            print(command.stderr)

            raise RuntimeError('Flask did not run successfully! :(')

    @staticmethod
    def get_file_num_lines(input_file):
        with open(input_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                pass
        return line_num

    def trim_flask_alm_output(self):
        """
        By default, Flask outputs the Alm file with a header.
        However, for this file to be read in by the C++ code we need to trim off this header.
        :return:
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

        :return: None
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

    def plot_flask_output(self):

        # Obtain the Planck colour-map used for plotting maps here
        planck_cmap = make_planck_colour_map()

        # Go through each redshift bin individually and make plots for each
        for z_bin in range(1, self.num_redshift_bins + 1):

            # Read in the generated shear map from Flask
            maps = hp.read_map(self.folder_path + 'Output-shear-map-f1z' + str(z_bin) + '.fits', verbose=False,
                               field=None)

            # Split the shae map into a convergence, shear1, and shear2 maps
            converg_map = maps[0]
            shear_map1 = maps[1]
            shear_map2 = maps[2]

            # converg_map = hp.read_map(self.folder_path + 'Output-shear-map-f1z' + str(z_bin) + '.fits', verbose=False)
            # shear_map1 = hp.read_map(self.folder_path + 'Output-shear-map-f1z' + str(z_bin) + '.fits', verbose=False, field=None)[1]
            # shear_map2 = hp.read_map(self.folder_path + 'Output-shear-map-f1z' + str(z_bin) + '.fits', verbose=False, field=None)[2]

            # Plot the convergence and shear1 maps
            hp.mollview(converg_map, cmap=planck_cmap, title='Convergence for redshift bin ' + str(z_bin))
            hp.graticule(verbose=False, alpha=0.6)
            hp.mollview(shear_map1, cmap=planck_cmap, title='Shear for redshift bin ' + str(z_bin))
            hp.graticule(verbose=False, alpha=0.6)

            # Use HealPy functions to turn the maps into Cl's
            start_time = time.time()
            converg_cl = hp.sphtfunc.anafast(converg_map, lmax=self.ell_max)
            print('This took {sec} seconds'.format(sec=time.time() - start_time))
            shear_cl1 = hp.sphtfunc.anafast(shear_map1, lmax=self.ell_max)
            shear_cl2 = hp.sphtfunc.anafast(shear_map2, lmax=self.ell_max)

            flask_cl_df = pd.read_csv(
                self.folder_path + 'Output-Reg-Cl-f1z' + str(z_bin) + 'f1z' + str(z_bin) + '.dat',
                header=None, names=['ell', 'Cl'], sep=r'\s+')

            # Plot various Cl's
            plt.figure(figsize=(12, 7))
            # plt.loglog(self.ells, self.ells * (self.ells + 1) * self.my_cls[z_bin - 1][2:] / (2 * np.pi), label='My Cl', color='yellow', lw=2.25)
            plt.loglog(self.ells, self.ells * (self.ells + 1) * converg_cl[2:] / (2 * np.pi), label=r'$C_\ell$ converg')
            plt.loglog(self.ells, self.ells * (self.ells + 1) * shear_cl1[2:] / (2 * np.pi), label=r'$C_\ell$ shear 1')
            plt.loglog(self.ells, self.ells * (self.ells + 1) * shear_cl2[2:] / (2 * np.pi), label=r'$C_\ell$ shear 2')
            plt.loglog(self.ells, self.ells * (self.ells + 1) * (shear_cl1[2:] + shear_cl2[2:]) / (2 * np.pi), label=r'$C_\ell$ shear 1 + 2', color='purple', lw=2)
            plt.loglog(self.ells, self.ells * (self.ells + 1) *
                       self.raw_c_ells['W' + str(z_bin) + 'x' + 'W' + str(z_bin)][2:] / (2 * np.pi),
                       label=r'$C_\ell$ input', lw=1.5, color='cyan')
            plt.loglog(flask_cl_df['ell'], flask_cl_df['ell'] * (flask_cl_df['ell'] + 1) * flask_cl_df['Cl'] / (2 * np.pi), lw=2, color='navy', label='Cl from Flask')
            plt.title(r'$C_\ell$ for redshift bin ' + str(z_bin))
            plt.xlabel(r'$\ell$')
            plt.ylabel(r'$\ell (\ell + 1) C_\ell / 2 \pi$')
            plt.legend()
            plt.tight_layout()

            # Plot various differences in the Cl's
            plt.figure(figsize=(12, 7))
            # plt.loglog(self.ells, np.abs(self.my_cls[z_bin - 1][2:] - converg_cl[2:]) / converg_cl[2:], label='Rel difference', color='blue', lw=1.5)
            plt.loglog(self.ells, np.abs(self.my_cls_cpp[z_bin - 1][2:] - converg_cl[2:]) / converg_cl[2:],
                       label='Rel. diff. cpp', color='purple', lw=1.5)
            plt.loglog(self.ells, np.abs(self.my_cls_cpp[z_bin - 1][2:] -
                                         self.raw_c_ells['W' + str(z_bin) + 'x' + 'W' + str(z_bin)][2:]) /
                       self.raw_c_ells['W' + str(z_bin) + 'x' + 'W' + str(z_bin)][2:],
                       label='Rel. diff. btw. input', color='blue', lw=1.5)
            plt.xlabel(r'$\ell$')
            plt.ylabel(r'$|C_{\ell}^{1} - C_{\ell}^{2} | / C_{\ell}^{2} $')
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

