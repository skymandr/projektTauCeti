#! /usr/bin/env python
"""
    surveyor.py
    - a python program, which shows an example use of projekt_anarres.py
    by letting the user survey a planet using different projections.

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

Usage:

    $ ./surveyor.py <rectangular map image>

Thie program lets you survey a map using the different projections implemented
in projekt_anarres.py. The followin projections are implemented:

    - azimuthal orthographic
    - azimuthal equal-area ("Lambert")
    - azimuthal equidistant
    - azimuthal stereographic
    - equirectangular

The map data is taken from an image file assumed by the program to contain
simple rectangular projection data. The program works best with png-format.

The azimuthal views only support views along the standard (equatorial)
parallel, while the rectangular view can be centered on any coordinate.

To update redraw meridians and parallels, press "Update". To reset the view,
press "Reset".

For more information on the different projections and options, please see
the documentation for projekt_anarres.py!


------
[0]:  http://www.gnu.org/licenses/gpl-3.0.html
"""


# Initial imports:
import os
import sys
import projekt_anarres as p
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, Cursor
from testgrids import make_testgrids as grids


class PlanetarySurveyor(object):
    """
    The PlanetarySurveyor creates a Matplotlib "widget" letting the user
    navigate map data loaded from an image file.
    """

    def __init__(self, filename):
        """
        Initialized with filename of image file containing the equirectangular
        map data.
        """

        self.filename = filename

        self.load_image()

        # Setup display:
        self.fig = plt.figure(1)
        self.ax = plt.subplot(111)
        plt.clf()
        plt.subplots_adjust(left=0.1, bottom=0.25)

        self.meridian = 90
        self.parallel = 90
        self.parallels = 16
        self.meridians = 16
        self.mode = "azimuthal"
        self.projection = "orthographic"

        self.setup_display()

        # Setup mouse interaction:
        self.click = self.display.figure.canvas.mpl_connect(
                        'button_press_event', self.mouseclick)

        self.cursor = Cursor(self.display.axes, useblit=True, color='red',
                             linewidth=1)

        # Setup axes:
        self.axes_step = plt.axes([0.13, 0.15, 0.60, 0.03])
        self.axes_meridians = plt.axes([0.13, 0.10, 0.60, 0.03])
        self.axes_parallels = plt.axes([0.13, 0.05, 0.60, 0.03])

        self.update_axes = plt.axes([0.79, 0.095, 0.08, 0.04])
        self.reset_axes = plt.axes([0.79, 0.05, 0.08, 0.04])

        self.radio_axes = plt.axes([0.88, 0.05, 0.11, 0.15])

        # Setup sliders:
        self.step = 22.5

        self.slider_step = Slider(self.axes_step, 'Step', 0, 90,
                                  valinit=self.step, valfmt='%2.1f')
        self.slider_meridians = Slider(self.axes_meridians, 'Meridians', 0,
                                       64, valinit=self.parallels,
                                       valfmt='%2d')
        self.slider_parallels = Slider(self.axes_parallels, 'Parallels', 0,
                                       64, valinit=self.parallels,
                                       valfmt='%2d')

        self.slider_step.on_changed(self.update)
        self.slider_meridians.on_changed(self.update)
        self.slider_parallels.on_changed(self.update)

        # Setup button(s):
        self.update_button = Button(self.update_axes, 'Update')
        self.update_button.on_clicked(self.update_display)

        self.reset_button = Button(self.reset_axes, 'Reset')
        self.reset_button.on_clicked(self.reset)

        # Setup radio buttons:
        self.radio = RadioButtons(self.radio_axes, ('ortho', 'eq.area',
                                  'eq.dist', 'stereo', 'rect'), active=0)
        self.radio.on_clicked(self.set_mode)
        self.projections = {"ortho": ("orthographic", "azimuthal"),
                            "eq.area": ("lambert", "azimuthal"),
                            "eq.dist": ("equidistant", "azimuthal"),
                            "stereo": ("stereographic", "azimuthal"),
                            "rect": ("rectangular", "rectangular")}

        # Almost ready:
        self.update()
        plt.show()

    def load_image(self):
        """
        Load and flatten specified image file. If this fails, the default map
        is loaded.  """

        try:
            map_image = plt.imread(self.filename)
        except IOError as e:
            print "Could not load file {0} ({1})".format(
                  self.filename, e.strerror)
            print "Using default image..."
            self.filename = "templates/nowwhat.png"
            map_image = plt.imread(self.filename)

        while len(map_image.shape) > 2:
            map_image = map_image.mean(-1)

        self.map_image = map_image

    def setup_display(self):
        """
        Setup parameters and map display.
        """

        self.R = 180
        self.padding = self.R / 10

        if self.mode == 'azimuthal':
            self.hemisphere = p.get_azimuthal_hemisphere(
                                self.map_image, self.meridian, 90,
                                self.R, self.projection, 0, self.padding)
            self.display = plt.imshow(self.hemisphere, cmap=plt.cm.gray,
                                      extent=[-1.5, 1.5, -1.5, 1.5])
            plt.axis([-1.5, 1.5, -1.5, 1.5])
            plt.axis('off')
        elif self.mode == 'rectangular':
            plt.axis([-180, 180, -90, 90])
            pass

        if self.meridians > 0 or self.parallels > 0:
            self.draw_graticules()

    def update_display(self, val=0):
        """
        Update map display.
        """

        if self.mode == 'azimuthal':
            self.hemisphere = p.get_azimuthal_hemisphere(
                                self.map_image, self.meridian, 90,
                                self.R, self.projection, 0, self.padding)
            ax = self.display.axes
            self.display.axes.cla()
            self.display = ax.imshow(self.hemisphere, cmap=plt.cm.gray,
                                      extent=[-1.5, 1.5, -1.5, 1.5])
            plt.axis([-1.5, 1.5, -1.5, 1.5])
            plt.axis('off')
        elif self.mode == 'rectangular':
            self.hemisphere = p.get_rectangular_hemisphere(
                                self.map_image, self.meridian, self.parallel,
                                False)
            ax = self.display.axes
            self.display.axes.cla()
            self.display = ax.imshow(self.hemisphere, cmap=plt.cm.gray,
                                     extent=[
                                     self.meridian - 90, self.meridian + 90,
                                     self.parallel - 90, self.parallel + 90])
            plt.axis([self.meridian - 90, self.meridian + 90,
                      self.parallel - 90, self.parallel + 90])

        self.fix_coordinates()

        if self.meridians > 0 or self.parallels > 0:
            self.draw_graticules()

        self.update()

    def update(self, val=0):
        """
        Update internal parameters from sliders, update coordiantes and draw.
        """

        if self.step != self.slider_step.val:
            self.step = np.round(self.slider_step.val / 0.5) * 0.5
            self.slider_step.set_val(self.step)

        if (self.meridians != self.slider_meridians.val or
                self.parallels != self.slider_parallels.val):

            self.meridians = np.round(self.slider_meridians.val
                                      ).astype(np.int)
            self.parallels = np.round(self.slider_parallels.val
                                      ).astype(np.int)

        self.fix_coordinates()

        plt.draw()

    def set_mode(self, val="ortho"):
        """
        Set projection and mode.
        """

        self.projection = self.projections[val][0]
        self.mode = self.projections[val][1]

        if self.mode == "azimuthal":
            self.parallel = 90

        self.update_display()

    def reset(self, event):
        """
        Reset widget
        """

        self.slider_step.reset()
        self.slider_meridians.reset()
        self.slider_parallels.reset()
        self.meridian = 90
        self.parallel = 90
        self.update()
        self.set_mode()

    def mouseclick(self, event):
        """
        Handle mouse navigation of map display for the different projections.
        """

        if event.inaxes == self.display.axes:
            if event.button == 1:
                if self.mode == "azimuthal":
                    self.meridian += self.step * np.round(event.xdata)
                elif self.mode == "rectangular":
                    self.parallel = 180 - np.round(event.ydata / 0.5) * 0.5
                    self.meridian = np.round(event.xdata / 0.5) * 0.5

                self.update()
                self.update_display()

    def fix_coordinates(self):
        """
        Fix coordinates so they comply to standard representation for maps.
        """

        self.parallel %= 180
        self.meridian %= 360

        self.display.axes.set_title("{0}: {1}".format(self.projection,
                                               self.get_coordinates()))

    def get_coordinates(self):
        """
        Return string representation of coordinates in N-S/E-W standard.
        """

        parallel = self.parallel
        meridian = self.meridian - 180

        if parallel > 90.0:
            parallel = 180.0 - parallel
        elif parallel < -90.0:
            parallel = -180.0 - parallel

        if meridian > 180.0:
            meridian = 360.0 - meridian
        elif meridian < -180.0:
            meridian = -360.0 - meridian

        if parallel >= 0:
            NS = "N"
        else:
            NS = "S"

        if meridian >= 0:
            EW = "E"
        else:
            EW = "W"

        return "{0} {1}, {2} {3}".format(np.abs(parallel), NS,
                                         np.abs(meridian), EW)

    def get_graticule(self):
        """
        Return resolution of the current graticule (distances between parallel
        and meridian lines).
        """

        try:
            dLat = 180.0 / self.parallels
        except ZeroDivisionError:
            dLat = None

        try:
            dLon = 360.0 / self.meridians
        except ZeroDivisionError:
            dLon = 0

        return dLat, dLon

    def draw_graticules(self):
        """
        Draw parallel and meridian lines using testgrids-module.
        """

        dLat, dLon = self.get_graticule()
        ax = self.display.axes

        if self.meridians > 0:
            if self.mode == "azimuthal":
                x_mer, z_mer = grids.get_meridians(self.meridian, dLon, 1,
                                                   1.5 / 1.1, self.projection)

                ax.plot(x_mer, z_mer, ':k', label="meridians")
                ax.axis('off')
            elif self.mode == "rectangular":
                mer_min = np.ceil((self.meridian - 90) / dLon) * dLon
                mer_max = np.floor((self.meridian + 90) / dLon) * dLon

                ax.set_xticks(np.arange(mer_min, mer_max + 1,  dLon))
                ax.grid(True)

        if self.parallels > 0:
            if self.mode == "azimuthal":
                x_par, z_par = grids.get_parallels(self.parallel, dLat, 1,
                                                   1.5 / 1.1, self.projection)

                ax.plot(x_par, z_par, ':k', label="parallels")
                ax.axis('off')
            elif self.mode == "rectangular":
                par_min = np.ceil((self.parallel - 90) / dLat) * dLat
                par_max = np.floor((self.parallel + 90) / dLat) * dLat

                ax.set_yticks(np.arange(par_min, par_max + 1, dLat))
                ax.grid(True)


def main():
    """
    Parse commandline arguments and start widget.
    """

    try:
        filename = sys.argv[1]
    except IndexError:
        filename = ""

    Surveyor = PlanetarySurveyor(filename)


if __name__ == "__main__":
    sys.exit(main())
