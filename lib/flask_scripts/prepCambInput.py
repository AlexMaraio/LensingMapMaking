#! /usr/bin/env python

"""
USAGE:   prepCambInput.py <CAMB_COV_FILE> <FIELDS_INFO_FILE> <FLASK_CL_PREFIX>
EXAMPLE: prepCambInput.py test08_scalCovCls.dat test08/fields-info.dat Cl-
OUTPUT:  <FLASK_CL_PREFIX>f1z1f1z1.dat, <FLASK_CL_PREFIX>f1z1f1z2.dat, ...

This script takes a CAMB angular power spectra output file with repeated columns and 
the FLASK FIELDS_INFO file <FIELDS_INFO_FILE> prepared by the camb2info.py script 
and writes a separate file for each of the power spectra in the <CAMB_COV_FILE> file 
containing two columns, ell and C(ell). Moreover, it removes the pi and
ell factors to return the pure C(ell)s.

Written by: Henrique S. Xavier, hsxavier@if.usp.br, 04/aug/2015.
"""

import numpy as np
import sys


# Function to get position of Cl in CAMB-sources Cov file:
def cov_position(f_1, f_2, n_fields):
    return (n_fields + 3) * (3 + f_1 - 1) + (3 + f_2)


if __name__ == '__main__':
    # Docstring output:
    if len(sys.argv) != 1 + 3:
        print(__doc__)
        sys.exit(1)

    # Get input:
    cambfile = sys.argv[1]
    infofile = sys.argv[2]
    outprefix = sys.argv[3]
    camb = np.loadtxt(cambfile, unpack=True)
    info = np.loadtxt(infofile)

    # Find out number of Fields:
    if len(info.shape) == 1:
        info = np.array([info])
    num_fields = info.shape[0]
    print('num_fields =', num_fields)

    # Prepare factors:
    ell = camb[0]
    fac = (2.0 * np.pi) / ell / (ell + 1)

    # Loop over fields x fields:
    for f1 in range(1, num_fields + 1):
        for f2 in range(f1, num_fields + 1):
            # Prepare Cls:
            Cl = np.transpose([ell, fac * camb[cov_position(f1, f2, num_fields)]])

            # Export:
            outfile = outprefix + 'f' + str(int(info[f1 - 1][0])) + 'z' + str(int(info[f1 - 1][1])) + 'f' + str(
                    int(info[f2 - 1][0])) + 'z' + str(int(info[f2 - 1][1])) + '.dat'
            print("Writing file " + outfile)
            np.savetxt(outfile, Cl, fmt=['%d', '%e'])

    print("Done.")


def split_files(input_file, output_file, num_redshift):
    print('Splitting the CAMB files')

    camb_output = np.loadtxt(input_file, unpack=True)

    ells = camb_output[0]

    # Loop over redshift x redshift bins:
    for z1 in range(1, num_redshift + 1):
        for z2 in range(z1, num_redshift + 1):
            # Obtain list of ells and c_ell values
            ell_c_l = np.transpose([ells, camb_output[cov_position(z1, z2, num_redshift)]])

            # Export:
            output_filename = output_file + 'f1' + 'z' + str(z1) + 'f1' + 'z' + str(z2) + '.dat'

            print("Writing file " + output_filename)
            np.savetxt(output_filename, ell_c_l, fmt=['%d', '%e'])

    print('Done outputting the files')
