#! /usr/bin/env python
import os
import sys
import projekt_anarres as p
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, Cursor
from testgrids import make_testgrids as grids


# TODO:
#   - grid
#   - add radio-buttons for:
#       * orthographic
#       * lambert
#       * equidistant
#       * stereographic
#       * rectangular
#       (the radio buttons should also set mode and parallel appropriately)
#   - implement those
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

        self.meridian = 90
        self.parallel = 90
        self.parallels = 0
        self.meridians = 0
        self.mode = "azimuthal"
        self.projection = "orthographic"

        self.setup_display()

        # Setup mouse interaction:
        self.click = self.display.figure.canvas.mpl_connect(
                        'button_press_event', self.mouseclick)

        self.cursor = Cursor(self.display.axes, useblit=True, color='red',
                             linewidth=1)

        # Setup axes:
        self.axes_step = plt.axes([0.15, 0.15, 0.60, 0.03])
        self.axes_meridians = plt.axes([0.15, 0.10, 0.60, 0.03])
        self.axes_parallels = plt.axes([0.15, 0.05, 0.60, 0.03])

        self.coord_axes = plt.axes([0.84, 0.14, 0.15, 0.04])
        self.update_axes = plt.axes([0.84, 0.095, 0.15, 0.04])
        self.reset_axes = plt.axes([0.84, 0.05, 0.15, 0.04])

        # Setup sliders:
        self.step = 22.5

        self.slider_step = Slider(self.axes_step, 'Step', 0, 90,
                                  valinit=self.step, valfmt='%2.1f')
        self.slider_meridians = Slider(self.axes_meridians, 'Meridians', 0,
                                       42, valinit=self.parallels,
                                       valfmt='%2d')
        self.slider_parallels = Slider(self.axes_parallels, 'Parallels', 0,
                                       42, valinit=self.parallels,
                                       valfmt='%2d')

        self.slider_step.on_changed(self.update)
        self.slider_meridians.on_changed(self.update)
        self.slider_parallels.on_changed(self.update)

        # Setup button(s):
        self.coords = Button(self.coord_axes, '{0}'.format(
                                              self.get_coordinates()))

        self.update_button = Button(self.update_axes, 'Update')
        self.update_button.on_clicked(self.update_display)

        self.reset_button = Button(self.reset_axes, 'Reset')
        self.reset_button.on_clicked(self.reset)

        plt.show()

    def load_image(self):
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
        if self.mode == 'azimuthal':
            self.hemisphere = p.get_azimuthal_hemisphere(
                                self.map_image, self.meridian, 90,
                                self.R, self.projection, 0, self.padding)
            ax = self.display.axes
            self.display = ax.imshow(self.hemisphere, cmap=plt.cm.gray,
                                      extent=[-1.5, 1.5, -1.5, 1.5])
            plt.axis([-1.5, 1.5, -1.5, 1.5])
            plt.axis('off')
        elif self.mode == 'rectangular':
            plt.axis([-180, 180, -90, 90])
            pass

        if self.meridians > 0 or self.parallels > 0:
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

        plt.draw()

    def reset(self, event):
        self.slider_step.reset()
        self.slider_meridians.reset()
        self.slider_parallels.reset()
        self.meridian = 90
        self.parallel = 90

    def mouseclick(self, event):
        if event.inaxes == self.display.axes:
            if event.button == 1:
                if self.mode == "azimuthal":
                    self.meridian += self.step * np.round(event.xdata)
                elif self.mode == "rectangular":
                    self.parallel = np.round(event.ydata / 0.5) * 0.5
                    self.meridian = np.round(event.xdata / 0.5) * 0.5

                self.fix_coordinates()

                self.update()
                self.update_display()

    def fix_coordinates(self):
            self.parallel %= 180
            self.meridian %= 360

            self.coords.label.set_text("{0}".format(self.get_coordinates()))

    def get_coordinates(self):
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
        try:
            dLat = np.round(180.0 / self.parallels).astype(np.int)
            dLon = np.round(360.0 / self.meridians).astype(np.int)

            return dLat, dLon
        except ZeroDivisionError:
            return None, None

    def draw_graticules(self):
        pass


def main():
    try:
        filename = sys.argv[1]
    except IndexError:
        filename = ""

    Surveyor = PlanetarySurveyor(filename)


if __name__ == "__main__":
    sys.exit(main())
