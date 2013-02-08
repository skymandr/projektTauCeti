#! /usr/bin/env python
"""
    surveyor_basemap.py
    - a reimplementation of surveyor.py, using the BaseMap package
    instead of projekt_anarres.py for drawing the projections.

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

    $ ./surveyor_basemap.py <rectangular map image>

Thie program lets you survey a map using some of the projections implemented
in projekt_anarres.py, but using BaseMap to calculate and draw the map. This
makes it more versetile than surveyor.py, of which it is a reimplementation,
but also slower. An aditional benefit is that this version works as a
standalone program, while surveyor.py requires the funtionality from
projekt_anarres.py and testgrids.  The following projections are available:

    - azimuthal orthographic
    - azimuthal equidistant
    - equirectangular

The map data is taken from an image file assumed by the program to contain
simple rectangular projection data. The program works best with png-format.

The views can be navigated along any parallel or meridian.

To update redraw meridians and parallels, press "Update". To reset the view,
press "Reset".

For more information on the different projections and options, please see
the documentation for projekt_anarres.py!


------
[0]:  http://www.gnu.org/licenses/gpl-3.0.html
"""


#Initial imports:
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, Cursor
from mpl_toolkits.basemap import Basemap


class PlanetarySurveyor(object):
    """
    The PlanetarySurveyor creates a Matplotlib "widget" letting the user
    navigate map data loaded from an image file.
    """

    def __init__(self, filename="templates/nowwhat.png"):
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
        plt.subplots_adjust(left=0.1, bottom=0.20)

        self.meridian = -90
        self.parallel = 22.5
        self.projection = "orthographic"
        self.parallels = 16
        self.meridians = 16

        self.setup_display()

        # Setup mouse interaction:
        self.click = self.hemisphere_axes.figure.canvas.mpl_connect(
                        'button_press_event', self.mouseclick)

        self.cursor = Cursor(self.hemisphere_axes, useblit=True, color='red',
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

        self.button = Button(self.reset_axes, 'Reset')
        self.button.on_clicked(self.reset)

        # Setup radio buttons:
        self.radio = RadioButtons(self.radio_axes, ('ortho', 'eq.dist',
                                                    'rect'), active=0)
        self.radio.on_clicked(self.set_mode)
        self.projections = {"ortho": ("orthographic", "ortho"),
                            "eq.dist": ("equidistant", "aeqd"),
                            "rect": ("rectangular", "cyl")}

        # Almost done:
        self.update()
        plt.show()

    def load_image(self):
        """
        Checks if the image file specified exists and is an image file. If this
        fails, the default map is loaded.
        """

        try:
            map_image = plt.imread(self.filename)
        except IOError as e:
            print "Could not load file {0} ({1})".format(
                  self.filename, e.strerror)
            print "Using default image..."
            self.filename = "templates/nowwhat.png"

    def setup_display(self):
        """
        Setup parameters and map display.
        """

        self.hemisphere_axes = plt.gca()
        self.R = 10000 * 360.0 / (2 * np.pi)
        self.width = 180 * 10000
        self.height = 180 * 10000

        if self.projection == 'orthographic':
            self.hemisphere = Basemap(projection="ortho", lon_0=self.meridian,
                              lat_0=self.parallel, resolution='c',
                              area_thresh=100000, ax=self.hemisphere_axes)
        elif self.projection == 'equidistant':
            self.hemisphere = Basemap(projection="aeqd", lon_0=self.meridian,
                              lat_0=self.parallel, resolution='c',
                              area_thresh=100000, rsphere=self.R,
                              ax=self.hemisphere_axes, width=self.width,
                              height=self.height)
        elif self.projection == 'rectangular':
            self.hemisphere = Basemap(projection="cyl", lon_0=0, lat_0=0,
                              resolution='c', area_thresh=100000,
                              ax=self.hemisphere_axes)

        self.hemisphere.warpimage(self.filename)
        self.hemisphere.drawmapboundary()
        self.draw_graticules()

    def update_display(self):
        """
        Update map display.
        """

        if self.projection == 'orthographic':
            self.hemisphere = Basemap(projection="ortho", lon_0=self.meridian,
                              lat_0=self.parallel, resolution='c',
                              area_thresh=100000, ax=self.hemisphere_axes)
        elif self.projection == 'equidistant':
            self.hemisphere = Basemap(projection="aeqd", lon_0=self.meridian,
                              lat_0=self.parallel, resolution='c',
                              area_thresh=100000, rsphere=self.R,
                              ax=self.hemisphere_axes, width=self.width,
                              height=self.height)
        elif self.projection == 'rectangular':
            self.hemisphere = Basemap(projection="cyl", lon_0=0, lat_0=0,
                              resolution='c', area_thresh=100000,
                              ax=self.hemisphere_axes)

        self.hemisphere_axes.cla()
        self.hemisphere.warpimage(self.filename)
        self.hemisphere.drawmapboundary()
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
        Set projection mode.
        """

        self.projection = self.projections[val][0]

        self.update_display()

    def reset(self, event):
        """
        Reset widget
        """

        self.slider_step.reset()
        self.slider_meridians.reset()
        self.slider_parallels.reset()
        self.meridian = 90
        self.parallel = 0
        self.update()

    def mouseclick(self, event):
        """
        Handle mouse navigation of map display for the different projections.
        """

        if event.inaxes == self.hemisphere_axes:
            if event.button == 1:
                if self.projection == "rectangular":
                    self.parallel = np.round(event.ydata / 0.5) * 0.5
                    self.meridian = np.round(event.xdata / 0.5) * 0.5
                else:
                    xlim = self.hemisphere_axes.get_xlim()
                    x = np.round(3 * (event.xdata -
                        0.5 * (xlim[1] - xlim[0])) / (xlim[1] - xlim[0]))

                    ylim = self.hemisphere_axes.get_ylim()
                    y = np.round(3 * (event.ydata -
                        0.5 * (ylim[1] - ylim[0])) / (ylim[1] - ylim[0]))

                    self.meridian += self.step * x
                    self.parallel += self.step * y

                    self.update_display()

                self.update()

    def fix_coordinates(self):
        """
        Fix coordinates so they comply to standard representation for maps.
        """

        if self.parallel > 90.0:
            self.parallel = 180.0 - self.parallel
        elif self.parallel < -90.0:
            self.parallel = -180.0 - self.parallel

        if self.meridian > 180.0:
            self.meridian = 360.0 - self.meridian
        elif self.meridian < -180.0:
            self.meridian = -360.0 - self.meridian

        self.hemisphere_axes.set_title("{0}: {1}".format(self.projection,
                                                  self.get_coordinates()))

    def get_coordinates(self):
        """
        Return string representation of coordinates in N-S/E-W standard.
        """

        parallel = np.abs(self.parallel)
        meridian = np.abs(self.meridian)

        if self.parallel >= 0:
            NS = "N"
        else:
            NS = "S"

        if self.meridian >= 0:
            EW = "E"
        else:
            EW = "W"

        return "{0} {1}, {2} {3}".format(parallel, NS, meridian, EW)

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
        Draw parallel and meridian lines.
        """

        parallel_step = 180.0 / self.parallels
        meridian_step = 360.0 / self.meridians

        if self.parallels > 0:
            min_parallel = -90
            max_parallel = 90 - parallel_step
            self.hemisphere.drawparallels(
                            np.arange(min_parallel, max_parallel + 1,
                            parallel_step), latmax=90 - parallel_step)

        if self.meridians > 0:
            min_meridian = -180
            max_meridian = 180 - meridian_step
            self.hemisphere.drawmeridians(
                            np.arange(min_meridian, max_meridian + 1,
                            meridian_step), latmax=90 - parallel_step)


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
