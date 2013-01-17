import projekt_anarres as p
import numpy as np
import matplotlib.pyplot as plt


def make_animation(the_map='templates/grid_double.png', saveas="animated",
                   frames=64, padding=10, padwith=0, R=256,
                   azikind='orthographic'):
    """
    HERE BE DOCSTRING!
    """
    
    for n in xrange(steps):
        hem = get_hemisphere(...)


def get_hemisphere(map_image, centre=90.0, R=256, azikind='orthographic',
                   padwith=0):
    """
    HERE BE DOCSTRING!
    """

    centre_coord = (np.round(
                   map_image.shape[1] * centre / 360.0))\
                   .astype(np.int)
    min_coord = (np.ceil(
                map_image.shape[1] * ((centre - 90.0) % 360.0) / 360.0))\
                .astype(np.int)
    max_coord = (np.ceil(
                map_image.shape[1] * ((centre - 90.0) % 360.0) / 360.0))\
                .astype(np.int)

    old_hemi, new_hemi = \
         p.really_make_azi_projection(rectangular, R, azikind, padwith)

    return new_hemi
    
