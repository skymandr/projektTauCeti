#! /usr/bin/env python
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, Cursor
from mpl_toolkits.basemap import Basemap


# TODO:
#   - make radiobuttons for the following projections:
#       * orthographic
#       * azimuthal equidistant
#       * equirectangular
#   - documentation


class PlanetarySurveyor(object):
    def __init__(self, filename):
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
        self.axes_step = plt.axes([0.15, 0.15, 0.60, 0.03])
        self.axes_meridians = plt.axes([0.15, 0.10, 0.60, 0.03])
        self.axes_parallels = plt.axes([0.15, 0.05, 0.60, 0.03])

        self.update_axes = plt.axes([0.84, 0.095, 0.15, 0.04])
        self.reset_axes = plt.axes([0.84, 0.05, 0.15, 0.04])

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

        # Almost done:
        self.update()
        plt.show()

    def load_image(self):
        try:
            map_image = plt.imread(self.filename)
        except IOError as e:
            print "Could not load file {0} ({1})".format(
                  self.filename, e.strerror)
            print "Using default image..."
            self.filename = "templates/nowwhat.png"

    def setup_display(self):
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
        if self.step != self.slider_step.val:
            self.step = np.round(self.slider_step.val / 0.5) * 0.5
            self.slider_step.set_val(self.step)

        if (self.meridians != self.slider_meridians.val or
                self.parallels != self.slider_parallels.val):

            self.meridians = np.int(self.slider_meridians.val)
            self.parallels = np.int(self.slider_parallels.val)

        self.fix_coordinates()

        plt.draw()

    def reset(self, event):
        self.slider_step.reset()
        self.slider_meridians.reset()
        self.slider_parallels.reset()
        self.meridian = 90
        self.parallel = 0
        self.update()

    def mouseclick(self, event):
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
        try:
            dLat = np.round(180.0 / self.parallels).astype(np.int)
            dLon = np.round(360.0 / self.meridians).astype(np.int)

            return dLat, dLon
        except ZeroDivisionError:
            return None, None

    def draw_graticules(self):
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
    try:
        filename = sys.argv[1]
    except IndexError:
        filename = ""

    Surveyor = PlanetarySurveyor(filename)


if __name__ == "__main__":
    sys.exit(main())
