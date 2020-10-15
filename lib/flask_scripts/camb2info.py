#! /usr/bin/env python

"""
USAGE:   camb2info.py <CAMB_INPUT> <FIELDS_INFO_FILENAME>
EXAMPLE: camb2info.py lcdm_10bins.ini lcdm_10bins/fields-info.dat

This script takes a CAMBsources input file and creates a file used by
FLASK describing the simulated fields and redshift bins (basically a
table of field and redshift bin IDs, mean values, shift parameters and
redshift ranges).

There are three methods for computing the shift parameter for convergence:
A table read from a file which is interpolated, the formula from
Hilbert, Hartlap & Schneider (2011), and a formula computed from FLASK
density line of sight integration. The latter is currently used (the other
ones are commented).

Written by Henrique S. Xavier, hsxavier@if.usp.br, 05/aug/2015.
"""


# Internal definitions:
GalMean = 0.0
KappaMean = 0.0
GalShift = 1.0
FixKappa = 0
KappaShift = 1.0
GalType = 1
KappaType = 2

# shiftfile  = "/home/skems/pkphotoz/prog/corrlnfields/data/Hilbert2011_shifts.dat"
# shiftfile  = "/home/skems/pkphotoz/prog/corrlnfields/data/k0_empty_LCDM_Om30.dat"

# Load shift file:
# ATTENTION!! Values not extrapolated: extremes are used!
# fileZ, fileSh = np.loadtxt(shiftfile, unpack=True)


# Convergence kappa shift formula from Hilbert, Hartlap & Schneider (2011)
def HilbertShift(z):
    return 0.008 * z + 0.029 * (z ** 2) - 0.0079 * (z ** 3) + 0.00065 * (z ** 4)


def XavierShift(z):
    a0 = 0.2
    s0 = 0.568591
    return a0 * (((z * s0) ** 2 + z * s0 + 1) / (z * s0 + 1) - 1)
