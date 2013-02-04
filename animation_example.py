"""
    animation_example.py:
    - a python module, which shows an example use of projekt_anarres.py

    (C) 2013 Andreas Skyman (skymandr@fripost.org)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see [0].


This module contains functions for gettiang azimuthal views centered on
arbitrary meridians from an equirectangular map, and for making and saving
frames for an animation based on such views

For more information on the different projections and options, please see
the documentation for projekt_anarres.py!


------
[0]:  http://www.gnu.org/licenses/gpl-3.0.html
"""


# Initial imports:

import projekt_anarres as p
import numpy as np
import matplotlib.pyplot as plt
import string


# Provided functions:

def make_animation(the_map='templates/grid_double.png', saveas="animated",
                   frames=64, padding=10, padwith=0, R=256,
                   azikind='orthographic'):
    """
    Fuction for creating frames in an animation of a rotating planet.
    """

    map_image = plt.imread(the_map)

    while len(map_image.shape) > 2:
        map_image = map_image.mean(-1)

    for n in xrange(frames):
        meridian = 90.0 + 360.0 * n / frames
        parallel = 90.0
        hemisphere = get_hemisphere(map_image, meridian, parallel,
                                    R, azikind, padwith)
        plt.imsave("{0}_{1}.png".format(saveas, string.zfill(n,
                   np.ceil(np.log10(frames)).astype(np.int))), hemisphere,
                   cmap=plt.cm.gray)


def get_hemisphere(map_image, meridian=90.0, parallel=90.0, R=256,
                   azikind='orthographic', padwith=0):
    """
    Function for getting an azimuthal view from a rectangular projection,
    centred on a particular meridian and parallel.

    (This function is also implemented in projekt_anarres.py.)
    """

    map_image = np.r_[map_image, np.flipud(
                np.c_[map_image[1:, map_image.shape[1] / 2:],
                      map_image[1:, : map_image.shape[1] / 2]])]

    meridian_coord = (np.round(
                   map_image.shape[1] * meridian / 360.0))\
                   .astype(np.int)
    min_meridian = (np.ceil(
                map_image.shape[1] * ((meridian - 90.0)) / 360.0))\
                .astype(np.int)
    max_meridian = (np.ceil(
                map_image.shape[1] * ((meridian + 90.0)) / 360.0))\
                .astype(np.int)

    parallel_coord = (np.round(
                   map_image.shape[0] * parallel / 360.0))\
                   .astype(np.int)
    min_parallel = (np.ceil(
                map_image.shape[0] * ((parallel - 90.0)) / 360.0))\
                .astype(np.int)
    max_parallel = (np.ceil(
                map_image.shape[0] * ((parallel + 90.0)) / 360.0))\
                .astype(np.int)

    Y, X = np.mgrid[min_parallel: max_parallel, min_meridian: max_meridian]

    X %= map_image.shape[1]
    Y %= map_image.shape[0]

    rectangular = map_image[(Y, X)]

    old_hemi, new_hemi = \
         p.really_make_azi_projection(rectangular, R, azikind, padwith)

    return new_hemi
