"""

Docstring again

"""

import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
import numpy as np
from shapely import geometry


from river import River
from ships import Ships

class UTM_Plotter:

    def __init__(self, river, ships, sim_set) -> None:
        
        print("Configuring UTM Plotter...")

        self.river = river
        self.ships = ships

       # Define colors:
        canvas_bg = "#bcb8b1"
        plot_bg = "#f4f3ee"

        # Get number of spawned ships on the plane
        #self.num_ships = ships.num_ships

        self.fig = plt.figure()
        self.fig.patch.set_facecolor(canvas_bg)

        gs = self.fig.add_gridspec(3,3)

        # Pane for the map
        self.ax1 = self.fig.add_subplot(gs[:, :-1])
        self.ax1.set_title('Some Map')

        # Make a contour plot with river coords as x and y, and stream velocity as z (why +10?)
        self.ax1.contourf(river.point_coords[:,:,0], river.point_coords[:,:,1],river.stream_vel + 10,
        cmap = cm.winter_r)
        self.ax1.set_facecolor(plot_bg)

        # Plot Vessels
        # xy coords are unpacked and plotted simultaneously for the utm transformed ship
        print("Placing vessels")
        for i in range(self.ships.num_ships):
            self.ax1.fill(*self.box_to_utm(self.ships.heading_box[i]).exterior.xy)
        
        # Panes for some metric to track 1
        self.ax2 = self.fig.add_subplot(gs[0, -1])
        self.ax2.set_title('Metric 1')
        self.ax2.set_facecolor(plot_bg)
        
        # Panes for some metric to track 2
        self.ax3 = self.fig.add_subplot(gs[1, -1])
        self.ax3.set_title('Metric 2')
        self.ax3.set_facecolor(plot_bg)
        
        # Panes for some metric to track 3
        self.ax4 = self.fig.add_subplot(gs[2, -1])
        self.ax4.set_title('Metric 3')
        self.ax4.set_facecolor(plot_bg)

    def render(self):
        plt.show()

    def box_to_utm(self, box):

        x,y = box.exterior.xy

        new_x = np.zeros(4)
        new_y = np.zeros(4)

        for i in range(4):
            new_x[i], new_y[i] = self.river.get_utm_position(x[i],y[i])

        new_polygon = geometry.Polygon([*zip(new_x,new_y)])

        return new_polygon


class Plotter:

    def __init__(self, river, ships, sim_set) -> None:
        

        print("Configuring Plotter...")

        self.river = river
        self.ships = ships

        # Define colors:
        canvas_bg = "#bcb8b1"
        plot_bg = "#f4f3ee"

        # Get number of spawned ships on the plane
        #self.num_ships = ships.num_ships

        self.fig = plt.figure()
        self.fig.patch.set_facecolor(canvas_bg)

        gs = self.fig.add_gridspec(4,1)

        # Pane for the map
        self.ax1 = self.fig.add_subplot(gs[0, :])
        self.ax1.set_title('Some Map')

        # TODO: Define a meshgrid for all x and y combinations in order to plot them
        # This is never used in the sim? (Viz.m:59)

        # Make a contour plot with river coords as x and y, and stream velocity as z (why +10?)
        print("Initialize Mesh...")
        xx,yy = np.meshgrid(np.arange(0,520,20),
        np.array([n*20+1 for n in range(len(self.river.point_coords))]))
        self.ax1.contourf(yy,xx,river.stream_vel + 10, cmap = cm.winter_r)
        self.ax1.set_facecolor(plot_bg)
        print("Done!")

        # Plot Vessels
        # xy coords are unpacked and plotted simultaneously for the utm transformed ship
        print("Placing vessels")
        for i in range(self.ships.num_ships):
            self.ax1.fill(*self.ships.heading_box[i].exterior.xy)
        
        # Panes for some metric to track 1
        self.ax2 = self.fig.add_subplot(gs[1, :])
        self.ax2.set_title('Metric 1')
        self.ax2.set_facecolor(plot_bg)
        
        # Panes for some metric to track 2
        self.ax3 = self.fig.add_subplot(gs[2, :])
        self.ax3.set_title('Metric 2')
        self.ax3.set_facecolor(plot_bg)
        
        # Panes for some metric to track 3
        self.ax4 = self.fig.add_subplot(gs[3,:])
        self.ax4.set_title('Metric 3')
        self.ax4.set_facecolor(plot_bg)

    def render(self):
        plt.show()



river = River()
ships = Ships(river=river, num_ships=5, ship_lengths= [100]*5, ship_widths= [10]*5,ship_mass=[1e5]*5,y_location=[150]*5)
viz = Plotter(river, ships, sim_set=None)
viz.render()
