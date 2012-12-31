"""
    projekt_anarres.py:
    - a python module for converting between azimuthal and
    equirectangular map projections

    (C) 2012 Andreas Skyman (skymandr@gmail.com)

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


This module contains functions for converting between orthographi, equidistant
and stereographic azimuthal [1, 2, 3], and equirectangular cylindrical [4] map
projections. While these two are arguably some of the cruder projections
available, they were chosen because they both make some sort of sense:

The orthographic projection is the simplest imaginable. At first, this appeared
to me to be the projection of choice in the "Hainish cycle" science fiction
works of Ursula K. Le Guin [5, 6], my interpetation of which was that this was
the standard for the space faring Ekumen: it made sense -- I argued -- since it
is the projection of the hemisphere onto the closest point on a a plane
bisecting the planet, corresponding to a view of the planet from very far away
(actually infinitely distant, but there you go...). Giving it a little more
thought, however, this projection makes little sense, because it discards so my
detail near the edge of the circle -- it would require a very special geography
of a planet for this to be a useful projection. The final desicion to move from
this kind of projection came when looking at maps, drawn by the author, of the
planet Gethen [7], where the polar regions were marked in such a way as to
indicate unambiguously, that the projection was not orthographic, but some more
advance azimuthal projection.

As an alternative to the orthographic projection, therefore, the module can
also handle equidistant and stereographic azimuthal projections. Those were a
little tricker implementing, but seem to yield much better results for all the
Hainish maps. The default projection is stereographic, but the qualitative
difference between the two is very slight.

The equirectangular projection is, perhaps, even less inspired [7], but this is
a standard often used by e.g. NASA in planetary survey maps [8]. It is the
simplest imaginable cylindrical projection, where all degrees are of draw with
the same length.


------
[0]: http://www.gnu.org/licenses/gpl-3.0.html
[1]: http://en.wikipedia.org/wiki/Orthographic_projection_(cartography)
[2]: http://en.wikipedia.org/wiki/Azimuthal_equidistant_projection
[3]: http://en.wikipedia.org/wiki/Stereographic_projection
[4]: http://en.wikipedia.org/wiki/Equirectangular_projection
[5]: http://www.ursulakleguin.com/FAQ.html#EkumenBooks
[6]: http://www.ursulakleguin.com/MenuContentsList.html#Illustrations
[7]: http://www.ursulakleguin.com/Maps/Map-Gethen.html
[8]: http://xkcd.com/977/
[9]: http://www.mapaplanet.org/
"""


# Initial imports:

import numpy as np
import matplotlib.pyplot as plt


# Provided functions:

def make_eqrec_projection(image='templates/stereographic.png',
                          azikind='stereographic',
                          latlong=(512, 512), cutoff=(0, 0)):
    """
    This function takes loads an image containing an azimuthal
    projection of (one hemisphere of a) a planetary body, and returns the
    corresponding equirectangular projection in a format specified by
    latlong.

    The azimuthal projection is assumed to be centered, and the padding
    around it to be small. The padding can be adjusted for by calling it with
    non-zero cutoff (unit is in pixels).

    The azimimuthal projection can be of either the "orthographic",
    "stereographic", or "equidistant" kind, as specified by the azikind
    parameter.

    If the image is not in grey scale, it is flattend first.

    Returns both old and new projection.
    """

    the_image = plt.imread(image)

    # Flatten image:
    while len(the_image.shape) > 2:
        the_image = the_image.mean(-1)

    # Prepare coordinates:
    Rz, Rx = the_image.shape[0] / 2, the_image.shape[1] / 2
    cox, coz = cutoff
    Lat, Long = (np.mgrid[0: 180: latlong[0] * 1j,
                          0: 180: latlong[1] * 1j]) * \
                 np.pi / 180.0

    # Make coordinate transform:
    X, Z = get_eqrec_coordinate_transform(Lat, Long, (Rx, Rz), (cox, coz),
                                          azikind)

    # Correct for nans:
    X = np.where(X < 2 * Rx, X, 0)
    X = np.where(X >= 0, X, 0)
    Z = np.where(Z < 2 * Rz, Z, 0)
    Z = np.where(Z >= 0, Z, 0)

    # Apply transform:
    new_image = the_image[Z, X]

    return the_image, new_image


def get_eqrec_coordinate_transform(Lat, Long, R, co, azikind='stereographic'):
    """
    Convenience function calculating the coordinate transform from
    azimuthal to equirectangular projection, where the azimuthal projection
    is either orthographic, stereographic (default), or equidistant.
    """

    Rx, Rz = R[0], R[1]
    cox, coz = co[0], co[1]

    if azikind.lower() == 'orthographic':
        X = (-np.round(np.sin(Long.T) * np.cos(Lat.T) \
                        * (Rx - 1 - cox)) + Rx).astype(np.int)
        Z = (-np.round(np.cos(Lat) * (Rz - 1 - coz)) + Rz).astype(np.int)
    elif azikind.lower() == 'equidistant':
        x = -np.sin(Long.T) * np.cos(Lat.T)
        z = -np.cos(Lat)
        y = np.sin(Long.T) * np.sin(Lat.T)

        d = np.arccos(y / 1)
        f = np.angle(x + 1j * z)

        X = (np.round(2 * (Rx - 1 - cox) * d * np.cos(f) \
            / np.pi + Rx)).astype(np.int)
        Z = (np.round(2 * (Rz - 1 - coz) * d * np.sin(f) \
            / np.pi + Rz)).astype(np.int)
    elif azikind.lower() == 'stereographic':
        x = -np.sin(Long.T) * np.cos(Lat.T)
        z = -np.cos(Lat)
        y = np.sin(Long.T) * np.sin(Lat.T)

        d = np.pi - np.arccos(y / 1)
        f = np.angle(x + 1j * z)

        rho = np.sin(d) / (1 - np.cos(d))
        the = f

        X = (np.round((Rx - 1 - cox) * rho * np.cos(the) \
            + Rx)).astype(np.int)
        Z = (np.round((Rz - 1 - coz) * rho * np.sin(the) \
            + Rz)).astype(np.int)
    else:
        try:
            raise OptionError(azikind)
        except:
            print 'Unknown option: {0} Assuming ortographic.'.format(azikind)
            return get_eqrec_coordinate_transform(Lat, Long, R, co,
                                                  'orthographic')

    return X, Z


def get_eqrec_map(east='templates/stereographic.png',
                  west='templates/stereographic.png',
                  azikind='stereographic', latlong=(720, 720), cutoff=(7, 1)):
    """
    This is a convenience function for creating a full equirectangular
    map from two hemisphere images, as described in make_eqrec_projection.
    """

    old_western, new_western = \
        make_eqrec_projection(west, azikind, latlong, cutoff)
    old_eastern, new_eastern = \
        make_eqrec_projection(east, azikind, latlong, cutoff)

    return np.hstack((new_western, new_eastern))


def draw_eqrec_map(the_map, dimensions=(-180, 180, -90, 90), origin='upper'):
    """
    This is a convenience function for drawing an equirectangular map
    with the right dimensions and and colour mode.

    NOTE: Different fileformats use different image coordinates for
    the bottom right pixel; some count the y-axis as starting in the
    bottom right corner, while others use top right corner as the origin.
    The script should work with png-files, but other formats (including jpg)
    may result in upside down images. If this happens, simply change origin
    from 'upper' to 'lower'.
    """

    plt.imshow(the_map, cmap=plt.cm.gray, extent=dimensions, origin=origin)

    plt.xticks(np.arange(dimensions[0], dimensions[1] + 1, 30))
    plt.yticks(np.arange(dimensions[2], dimensions[3] + 1, 30))


def make_azi_projection(image='templates/grid_double.png',
                        both_hemispheres=True, R=256,
                        azikind='stereographic', padwith=0):
    """
    This is a wrapper function for really_make_azi_projection, to make
    it simpler to automatically make two hemispheres in one go from a single
    source image.

    The function loads a specified image and flattens it if it is not in
    grey scale already. It then calls recall_make_azi_projection the
    appropriate number of times to perform the azimuthal projection. The
    azimuthal porjection can be of either the orthographic, the sterographic
    (default) or the equidistant kind.

    The parameter padwith decides the hue of the area outside the map.
    It should be between 0 (black) and 1 (white).

    Returns either one or two hemispheres, depending of the value of
    both_hemispheres.
    """

    the_map = plt.imread(image)

    # Flatten image:
    while len(the_map.shape) > 2:
        the_map = the_map.mean(-1)

    if both_hemispheres:
        shape = the_map.shape
        old_west, new_west = \
            really_make_azi_projection(the_map[:, : shape[1] / 2], R,
                                       azikind, padwith)
        old_east, new_east = \
            really_make_azi_projection(the_map[:,   shape[1] / 2:], R,
                                       azikind, padwith)
        return new_west, new_east
    else:
        old_hemi, new_hemi = really_make_azi_projection(the_map, R)
        return new_hemi


def really_make_azi_projection(hemisphere, R=256,
                               azikind='stereographic', padwith=0):
    """
    From supplied image data, assumed to be one hemisphere of an
    equirectangular projection, transforms this into an azimuthal
    projection. The resulting array has the shape (2*R+1, 2*R+1).

    Can handle both orthographic, stereographic (default), and
    equirectangular azimuthal projections, as specified by the
    parameter azikind.

    The parameter padwith decides the hue of the area outside the map.
    It should be between 0 (black) and 1 (white).

    Returns both orignal hemisphere data and azimuthal projection.
    """

    S = np.min(hemisphere.shape)

    # Prepare coordinates:
    Z, X = np.mgrid[-1.0: 1.0: R * 1j, -1.0: 1.0: R * 1j]
    Y = -np.sqrt(1.0 - Z ** 2 - X ** 2)

    # Make coordinate transform:
    Lat, Long = get_azi_coordinate_transform(X, Y, Z, S, azikind)

    # Correct for nans:
    Lat = np.where(Lat <= S, Lat, 0)
    Lat = np.where(Lat >= 0, Lat, 0)
    Long = np.where(Long <= S, Long, 0)
    Long = np.where(Long >= 0, Long, 0)

    # Apply transform:
    azimuthal = hemisphere[(Lat, Long)]
    azimuthal = np.where(X ** 2 + Y ** 2 <= 1.0, azimuthal, padwith)

    return hemisphere, azimuthal


def get_azi_coordinate_transform(X, Y, Z, S, azikind='stereographic'):
    """
    Convenience function calculating the coordinate transform from
    equirectangular to azimuthal projection, where the azimuthal projection
    is either orthographic, stereographic (default), or equidistant.
    """

    if azikind.lower() == 'orthographic':
        Lat  = (np.round(np.arccos(-Z / 1) * (S - 1) / np.pi)).astype(np.int)
        Long = (np.round(((np.arctan(Y / X) + np.pi) % np.pi) \
                          * (S - 1) / np.pi)).astype(np.int)
    elif azikind.lower() == 'equidistant':
        d = np.sqrt(X ** 2 + Z ** 2) * np.pi / 2
        f = np.angle(X + 1j * Z)

        x = -np.cos(d)
        y = np.sin(d) * np.cos(f)
        z = -np.sin(d) * np.sin(f)

        Lat  = (np.round(np.arccos(z / 1) * \
                         (S - 1) / np.pi))\
                         .astype(np.int)
        Long = (np.round((np.arctan(x / y) % np.pi) * \
                         (S - 1) / np.pi))\
                         .astype(np.int)
    elif azikind.lower() == 'stereographic':
        d = 2 * np.arctan(1 / np.sqrt(X ** 2 + Z ** 2))
        f = np.angle(X + 1j * Z)
        d -= np.pi

        x = np.cos(d)
        y = np.sin(d) * np.cos(f)
        z = np.sin(d) * np.sin(f)

        Lat  = (np.round(np.arccos(z / 1) * \
                         (S - 1) / np.pi))\
                         .astype(np.int)
        Long = (np.round((np.arctan(x / y) % np.pi) * \
                         (S - 1) / np.pi))\
                         .astype(np.int)
    else:
        try:
            raise OptionError(azikind)
        except:
            print 'Unknown option: {0} Assuming stereographic.'.format(azikind)
            return get_azi_coordinate_transform(X, Y, Z, S, 'stereographic')

    return Lat, Long


def get_azi_map(the_map='templates/grid_double.png', padding=10, padwith=0,
                R=256, azikind='stereographic'):
    """
    This is a convenience function for creating a composite map of both
    hemispheres with azimuthal projection from a single image, as described in
    *make_azi_projection. The default azimuthal projection is stereographic,
    with orthoographic and equidistant as alternatives.

    Padding is added around the azimutahl hemispheres with the amount
    requested. The colour of the padding is black (padwith=0) by default, but
    can be any grey scale value between 0 and 1.

    Returns composite azimuthal map.
    """

    # Get projection data:
    western, eastern = make_azi_projection(the_map, True, R, azikind, padwith)

    # Apply padding:
    new_western = np.zeros((western.shape[0] + 2 * padding,
                            western.shape[1] + 2 * padding)) + padwith
    new_eastern = np.zeros((eastern.shape[0] + 2 * padding,
                            eastern.shape[1] + 2 * padding)) + padwith

    new_western[padding: -padding, padding:-padding] = western
    new_eastern[padding: -padding, padding:-padding] = eastern

    return np.hstack((new_western, new_eastern))


def draw_azi_map(the_map, origin='upper'):
    """
    This is a convenience function for drawing an azimuthal map
    with the right dimensions and and colour mode.

    NOTE: Different fileformats use different image coordinates for
    the bottom right pixel; some count the y-axis as starting in the
    bottom right corner, while others use top right corner as the origin.
    The script should work with png-files, but other formats (including jpg)
    may result in upside down images. If this happens, simply change origin
    from 'upper' to 'lower'.
    """

    plt.imshow(the_map, cmap=plt.cm.gray, origin=origin)

    plt.axis('off')


def test_projection(the_image='templates/stereographic.png',
                    azikind='stereographic',
                    latlong=(512, 512), cutoff=(0, 0), origin='upper'):
    """
    This is a function for testing the efficacy of the forward and reverse
    transforms.
    """

    im, nim = make_eqrec_projection(the_image, azikind, latlong, cutoff)

    if im.shape[0] != im.shape[1]:
        print "'Diff' will only work properly for square maps!"
        print "Will cheat now..."
        the_shape = np.array(im.shape)
        R = the_shape.min()
        hem, azi = really_make_azi_projection(nim, R, azikind)
        im = im[:R, :R]
    else:
        R = im.shape[0]
        hem, azi = really_make_azi_projection(nim, R, azikind)

    plt.subplot(221)
    plt.title('original:')
    plt.imshow(im, cmap=plt.cm.gray, origin=origin)

    plt.subplot(222)
    plt.title('{0}'.format(azikind))
    plt.imshow(nim, cmap=plt.cm.gray, origin=origin)

    plt.subplot(223)
    plt.title('reverse {0}:'.format(azikind))
    plt.imshow(azi, cmap=plt.cm.gray, origin=origin)

    plt.subplot(224)
    plt.title('diff:')
    plt.imshow(azi - im, cmap=plt.cm.gray, origin=origin)


def get_rectangular_grid(spacing=30, size=720, both=False):
    """
    Help function for creating simple rectangular grid, for debugging
    and/or testing the different projections.

    Returns grid array.
    """

    if both:
        M, P = np.mgrid[0: size + 1: spacing, 0: 2 * size + 1: spacing]
        grid = np.ones((size + 1, 2 * size + 1))
    else:
        M, P = np.mgrid[0: size + 1: spacing, 0: size + 1: spacing]
        grid = np.ones((size + 1, size + 1))

    grid[M, :] = 0
    grid[:, P] = 0

    return grid


class OptionError(Exception):
    """
    Error raised if invalid option is passed to a function.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
