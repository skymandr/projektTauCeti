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
#   - file assignment on startup
#   - documentation


class PlanetarySurveyor(object):
    def __init__(self, filename):
        self.filename = filename

        map_image = plt.imread(filename)
        while len(map_image.shape) > 2:
            map_image = map_image.mean(-1)

        self.map_image = map_image

        # Setup display:
        self.fig = plt.figure(1)
        self.ax = plt.subplot(111)
        plt.clf()
        plt.subplots_adjust(left=0.1, bottom=0.20)

        self.meridian = 90
        self.parallel = 90
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
        self.reset_axes = plt.axes([0.84, 0.05, 0.15, 0.04])
        self.coord_axes = plt.axes([0.84, 0.14, 0.15, 0.04])

        # Setup sliders:
        self.step = 22.5
        self.parallels = 0
        self.meridians = 0

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
        self.button = Button(self.reset_axes, 'Reset')
        self.button.on_clicked(self.reset)

        self.coords = Button(self.coord_axes, '{0}'.format(
                                              self.get_coordinates()))

        plt.show()

    def setup_display(self):
        self.R = 180
        self.padding = self.R / 10
        print self.padding
        if self.mode == 'azimuthal':
            self.hemisphere = p.get_azimuthal_hemisphere(
                                self.map_image, self.meridian, 0,
                                self.R, self.projection, 0, self.padding)
            self.display = plt.imshow(self.hemisphere, cmap=plt.cm.gray,
                                      extent=[-1.5, 1.5, -1.5, 1.5])
            plt.axis([-1.5, 1.5, -1.5, 1.5])
            plt.axis('off')
        elif self.mode == 'rectangular':
            plt.axis([-180, 180, -90, 90])
            pass

    def update_display(self):
        if self.mode == 'azimuthal':
            self.hemisphere = p.get_azimuthal_hemisphere(
                                self.map_image, self.meridian, 0,
                                self.R, self.projection, 0, self.padding)
            ax = self.display.axes
            self.display = ax.imshow(self.hemisphere, cmap=plt.cm.gray,
                                      extent=[-1.5, 1.5, -1.5, 1.5])
            plt.axis([-1.5, 1.5, -1.5, 1.5])
            plt.axis('off')
        elif self.mode == 'rectangular':
            plt.axis([-180, 180, -90, 90])
            pass

    def update(self, val):
        if self.step != self.slider_step.val:
            self.step = np.round(self.slider_step.val / 0.5) * 0.5
            self.slider_step.set_val(self.step)

        if (self.meridians != self.slider_meridians.val or
                self.parallels != self.slider_parallels.val):

            self.meridians = np.int(self.slider_meridians.val)
            self.parallels = np.int(self.slider_parallels.val)

            if self.meridians > 0 or self.parallels > 0:
                self.draw_graticules()
            else:
                self.update_display()

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

                self.update(0)
                self.update_display()

    def fix_coordinates(self):
            if self.parallel > 90.0:
                self.parallel = 180.0 - self.parallel
            elif self.parallel < -90.0:
                self.parallel = -180.0 - self.parallel

            if self.meridian > 180.0:
                self.meridian = 360.0 - self.meridian
            elif self.meridian < -180.0:
                self.meridian = -360.0 - self.meridian

            self.coords.label.set_text("{0}".format(self.get_coordinates()))

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
        self.update_display()


def main():
    filename = 'new/anarres_small_eq_2_grey.png'
    print filename
    Surveyor = PlanetarySurveyor(filename)


if __name__ == "__main__":
    sys.exit(main())
