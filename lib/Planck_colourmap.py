"""

This file allows for the use of the Planck colour-map style for use in HealPy maps.

"""


from matplotlib.colors import ListedColormap
import numpy as np


def make_planck_colour_map():
    """
    Function to return the Planck-style colour map

    For example, this can be used as as colormap directly with hp.mollview(map, cmap=make_planck_colour_map())

    Note: Requires the correct location of the Planck_Parchment_RGB.txt file

    :return: Matplotlib colour map class
    """

    # Loads in data and imports to colour map
    cmap = ListedColormap(np.loadtxt("./resources/Planck_Parchment_RGB.txt")/255.)

    # Colour of missing pixels
    cmap.set_bad("gray")

    # Sets the colour of background, if necessary
    cmap.set_under("white")

    return cmap
