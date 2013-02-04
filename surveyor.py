import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
import projekt_anarres as p


class PlanetarySurveyor(object):
    def __init__(self, filename):
        self.filename = filename

        map_image = plt.imread(filename)
        while len(map_image.shape) > 2:
            map_image = map_image.mean(-1)

        self.map_image = map_image

        self.fig = plt.figure(1)
        self.ax = plt.subplot(111)
        plt.clf()
        plt.subplots_adjust(left=0.1, bottom=0.20)
        self.meridian = 90
        self.parallel = 90
        self.R = 360
        self.padding = 36
        hemisphere = p.get_hemisphere(self.map_image, self.meridian,
                                      self.parallel, self.R,
                                      padding=self.padding)
        self.display = plt.imshow(hemisphere, cmap=plt.cm.gray,
                                  extent=[-90, 90, -90, 90])
        plt.axis([-90, 90, -90, 90])
        plt.axis('off')
        self.axes_step = plt.axes([0.15, 0.15, 0.60, 0.03])
        self.axes_meridians = plt.axes([0.15, 0.10, 0.60, 0.03])
        self.axes_parallels = plt.axes([0.15, 0.05, 0.60, 0.03])

        self.step = 90
        self.parallels = 0
        self.meridians = 0

        self.slider_step = Slider(self.axes_step, 'Step', 0, 90,
                                  valinit=self.step, valfmt='%2d')
        self.slider_meridians = Slider(self.axes_meridians, 'Meridians', 0,
                                       360, valinit=self.parallels,
                                       valfmt='%2d')
        self.slider_parallels = Slider(self.axes_parallels, 'Parallels', 0,
                                       180, valinit=self.parallels,
                                       valfmt='%2d')

        self.slider_step.on_changed(self.update)
        self.slider_meridians.on_changed(self.update)
        self.slider_parallels.on_changed(self.update)

        self.resetax = plt.axes([0.8, 0.025, 0.1, 0.04])
        self.button = Button(self.resetax, 'Reset')
        self.button.on_clicked(self.reset)

        plt.show()

    def update(self, val):
        self.step = np.int(self.slider_step.val)
        self.meridians = np.int(self.slider_meridians.val)
        self.parallels = np.int(self.slider_parallels.val)
        self.hemisphere = p.get_hemisphere(self.map_image, self.meridians,
                                           self.parallels, self.R,
                                           padding=self.padding)
        self.display.set_data(self.hemisphere)
        #print self.get_coordinates()
        #print self.get_graticule()
        plt.draw()

    def reset(self, event):
        self.slider_step.reset()
        self.slider_meridians.reset()
        self.slider_parallels.reset()

    def get_coordinates(self):
        meridian = self.meridians % 360
        parallel = self.parallels % 360

        return parallel, meridian

    def get_graticule(self):
        try:
            dLat = np.round(180.0 / self.parallels).astype(np.int)
            dLon = np.round(360.0 / self.meridians).astype(np.int)

            return dLat, dLon
        except ZeroDivisionError:
            return None, None

    def draw_meridians(self):
        pass

    def draw_parallels(self):
        pass

filename = 'new/anarres_small_eq_2_grey.png'

Surveyor = PlanetarySurveyor(filename)
