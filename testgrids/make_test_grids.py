"""
    make-test_grids.py
    - a python module for testing and showcasing the projections
    implemented in projekt_anarres.py

    (C) 2012 Andreas Skyman (skymandr@fripost.org)

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


This module contains functions for converting between orthographic, equidistant
and stereographic azimuthal [1, 2, 3], and equirectangular cylindrical [4] map
projections, with the aim of showing the difference between them.

Please see projekt_anarres.py for more detailed information!


------
[0]: http://www.gnu.org/licenses/gpl-3.0.html
[1]: http://en.wikipedia.org/wiki/Orthographic_projection_(cartography)
[2]: http://en.wikipedia.org/wiki/Azimuthal_equidistant_projection
[3]: http://en.wikipedia.org/wiki/Stereographic_projection
[4]: http://en.wikipedia.org/wiki/Equirectangular_projection
"""

import numpy as np
import matplotlib.pyplot as plt
from  matplotlib import rc
import projekt_anarres as p


def get_eqrec_coordinate_transform(Lat, Lon, R=1.0, azikind='stereographic'):
    """
    Convenience function calculating the coordinate transform from
    azimuthal to equirectangular projection, where the azimuthal projection
    is either orthographic, stereographic (default), or equidistant.
    """

    if azikind.lower() == 'orthographic':
        X = - R * np.sin(Lon.T) * np.cos(Lat.T)
        Z = - R * np.cos(Lon.T)
    elif azikind.lower() == 'equidistant':
        x = -np.sin(Lon.T) * np.cos(Lat.T)
        z = -np.cos(Lon.T)
        y = np.sin(Lon.T) * np.sin(Lat.T)

        d = np.arccos(y / 1)
        f = np.angle(x + 1j * z)

        X = R * 2 * d * np.cos(f) / np.pi
        Z = R * 2 * d * np.sin(f) / np.pi
    elif azikind.lower() == 'stereographic':
        x = -np.sin(Lon.T) * np.cos(Lat.T)
        z = -np.cos(Lon.T)
        y = np.sin(Lon.T) * np.sin(Lat.T)

        d = np.pi - np.arccos(y / 1)
        f = np.angle(x + 1j * z)

        rho = np.sin(d) / (1 - np.cos(d))
        the = f

        X = R * rho * np.cos(the)
        Z = R * rho * np.sin(the)
    else:
        try:
            raise OptionError(azikind)
        except:
            print 'Unknown option: {0} Assuming ortographic.'.format(azikind)
            return get_eqrec_coordinate_transform(Lat, Lon, R, 'orthographic')

    return X, Z


def get_azi_coordinate_transform(X, Z, S, azikind='stereographic'):
    """
    Convenience function calculating the coordinate transform from
    equirectangular to azimuthal projection, where the azimuthal projection
    is either orthographic, stereographic (default), or equidistant.
    """

    Y = -np.sqrt(1.0 - X ** 2 - Z ** 2)

    if azikind.lower() == 'orthographic':
        Lat = np.arccos(-Z / 1) * S / np.pi
        Lon = (np.arctan(Y / X) + np.pi) % np.pi * S / np.pi
    elif azikind.lower() == 'equidistant':
        d = np.sqrt(X ** 2 + Z ** 2) * np.pi / 2
        f = np.angle(X + 1j * Z)

        x = -np.cos(d)
        y = np.sin(d) * np.cos(f)
        z = -np.sin(d) * np.sin(f)

        Lat = np.arccos(z / 1) * S / np.pi
        Lon = np.arctan(x / y) % np.pi * S / np.pi
    elif azikind.lower() == 'stereographic':
        d = 2 * np.arctan(1 / np.sqrt(X ** 2 + Z ** 2))
        f = np.angle(X + 1j * Z)
        d -= np.pi

        x = np.cos(d)
        y = np.sin(d) * np.cos(f)
        z = np.sin(d) * np.sin(f)

        Lat  np.arccos(z / 1) * S / np.pi
        Lon = np.arctan(x / y) % np.pi * S / np.pi
    else:
        try:
            raise OptionError(azikind)
        except:
            print 'Unknown option: {0} Assuming stereographic.'.format(azikind)
            return get_azi_coordinate_transform(X, Y, Z, S, 'stereographic')

    return Lat, Lon


def make_test(graticule=15, resolution=1, fign=1, clearit=True,
              projections=('orthographic', 'stereographic', 'equidistant')):
    """
    Function for producing a figure showcasing the differences between
    the selected projections.
    """

    fig = get_new_figure(fign)
    if clearit:
        fig.clf()

    Lat1, Lon1 = np.mgrid[0: 181: resolution, 0: 181: graticule] * np.pi / 180
    Lat2, Lon2 = np.mgrid[0: 181: graticule, 0: 181: resolution] * np.pi / 180

    N = len(projections)
    mc = 1 - 1e-6

    for s, src in enumerate(projections):
        x1, z1 = get_eqrec_coordinate_transform(Lat1, Lon1, 1.0, src)
        x2, z2 = get_eqrec_coordinate_transform(Lat2, Lon2, 1.0, src)
        plt.subplot(N, N + 1, (N + 1) * s + 1)
        plt.plot(x1.T, z1.T, 'r', label='Parallels')
        plt.plot(x2, z2, 'k', label='Meridians')

        ax = plt.gca()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        plt.axis('equal')
        plt.axis([-1.04, 1.04, -1.04, 1.04])

        plt.ylabel(src)
        if s == 0:
            plt.title("original:")

        for d, dest in enumerate(projections):
            la1, lo1 = get_azi_coordinate_transform(x1 * mc, z1 * mc,
                                                    180, dest)
            la2, lo2 = get_azi_coordinate_transform(x2 * mc, z2 * mc,
                                                    180, dest)

            plt.subplot(N, N + 1, (N + 1) * s + 2 + d)
            plt.plot(lo1.T, la1.T, 'r', label='Parallels')
            plt.plot(lo2, la2, 'k', label='Meridians')

            ax = plt.gca()
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            plt.axis('equal')
            plt.axis([-4, 184, -4, 184])

            if s == 0:
                plt.title("assumption:\n" + dest, size='medium')


def get_new_figure(fign):
    rc('text', usetex=True)
    rc('font', family='serif')
    rc('font', size=8.0)
    rc('lines', markersize=5)
    rc('lines', linewidth=0.6)
    rc('lines', markeredgewidth=1.0)
    rc('legend', fontsize='x-small')
    rc('legend', columnspacing=1.0)
    rc('savefig', dpi=600)
    f = plt.figure(fign, figsize=(8, 5), dpi=150)
    f.subplots_adjust(left=0.15, bottom=0.15)
    return f


def old_test_grids():
    """
    Function for making a comparison of the different projections in
    bit map from. Reduntant.
    """

    things = ['orthographic', 'stereographic', 'equidistant']
    for s, src in enumerate(things):
        for d, dest in enumerate(things):
            im, nim = p.make_eqrec_projection('tempaltes/{0}.png'.format(src),
                                              azikind=dest, cutoff=(0, 0))
            plt.imsave('testgrids/{0}_{1}'.format(src, dest), nim,
                       cmap=cm.gray)
            plt.subplot(3, 3, 3 * s + d + 1)
            plt.imshow(nim, cmap=cm.gray)
            if d == 0:
                plt.ylabel("original:\n" + src)
            if s == 0:
                plt.title("assumption:\n" + dest, size='small')
            ax = plt.gca()
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xticklabels([])
            ax.set_yticklabels([])


class OptionError(Exception):
    """
    Error raised if invalid option is passed to a function.
    """

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
