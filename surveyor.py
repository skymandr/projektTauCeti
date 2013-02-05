#! /usr/bin/env python
import os
import sys
import projekt_anarres as p
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons, Cursor
from testgrids import make_testgrids as grids

# TODO:
#   - think about implementing arbitrary azimuthal
#   - rectangular mode
#   - radio buttons for picking mode
#   - show grid
#   - fix map mirror bug in projekt_anarres


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

        self.setup_display()

        # Setip mouse interaction:
        self.click = self.display.figure.canvas.mpl_connect(
                        'button_press_event', self.mouseclick)

        self.cursor = Cursor(self.display.axes, useblit=True, color='red',
                             linewidth=1)

        # Setup axes:
        self.axes_step = plt.axes([0.15, 0.15, 0.60, 0.03])
        self.axes_meridians = plt.axes([0.15, 0.10, 0.60, 0.03])
        self.axes_parallels = plt.axes([0.15, 0.05, 0.60, 0.03])
        self.reset_axes = plt.axes([0.82, 0.05, 0.15, 0.04])
        self.coord_axes = plt.axes([0.82, 0.14, 0.15, 0.04])

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
        self.R = 360
        self.padding = self.R / 10
        if self.mode == 'azimuthal':
            self.hemisphere = p.get_azimuthal_hemisphere(
                                self.map_image, self.meridian, self.parallel,
                                self.R, padding=self.padding)
            self.display = plt.imshow(self.hemisphere, cmap=plt.cm.gray,
                                      extent=[-1.5, 1.5, -1.5, 1.5])
            plt.axis([-1.5, 1.5, -1.5, 1.5])
            plt.axis('off')
        elif self.mode == 'rectangular':
            pass

    def update_display(self):
        if self.mode == 'azimuthal':
            self.hemisphere = p.get_azimuthal_hemisphere(
                                self.map_image, self.meridian, self.parallel,
                                self.R, padding=self.padding)
            ax = self.display.axes
            self.display = ax.imshow(self.hemisphere, cmap=plt.cm.gray,
                                      extent=[-1.5, 1.5, -1.5, 1.5])
            plt.axis([-1.5, 1.5, -1.5, 1.5])
            plt.axis('off')
        elif self.mode == 'rectangular':
            pass

    def update(self, val):
        if self.step != self.slider_step.val:
            self.step = np.round(self.slider_step.val / 0.5) * 0.5
            self.slider_step.set_val(self.step)
        self.meridians = np.int(self.slider_meridians.val)
        self.parallels = np.int(self.slider_parallels.val)
        
        self.update_display()

        if self.meridians > 0 or self.parallels > 0:
            self.draw_graticules()

        #print self.get_coordinates()
        #print self.get_graticule()

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
                x, y = np.round(event.xdata), np.round(event.ydata)
                if self.mode == "azimuthal":
                    self.meridian += self.step * x
                    #self.parallel -= self.step * y
                elif self.mode == "rectangular":
                    self.meridian = x
                    self.parallel = -y

                self.meridian %= 360.0
                self.parallel %= 360.0
                self.coords.label.set_text("{0}".format(
                                           self.get_coordinates()))
                self.update(0)

    def get_coordinates(self):
        parallel = self.parallel
        meridian = self.meridian

        return parallel, meridian

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
    filename = 'new/anarres_small_eq_2_grey.png'
    print filename
    Surveyor = PlanetarySurveyor(filename)


if __name__ == "__main__":
    sys.exit(main())
