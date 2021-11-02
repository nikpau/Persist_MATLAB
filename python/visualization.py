"""

Docstring again

"""

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
from river import River

class Plotter:

    def __init__(self, river, ships, sim_set) -> None:
            
        # Get number of spawned ships on the plane
#        self.num_ships = ships.num_ships

        self.fig = plt.figure()
 
        gs = self.fig.add_gridspec(3,3)

        # Pane for the map
        self.ax1 = self.fig.add_subplot(gs[:, :-1])
        self.ax1.set_title('Some Map')

        # TODO: Define a meshgrid for all x and y combinations in order to plot them
        # This is never used in the sim? (Viz.m:59)

        # Make a contour plot with river coords as x and y, and stream velocity as z (why +10?)
        self.ax1.contour(river.point_coords[:,:,1], river.point_coords[:,:,0],river.stream_vel + 10)
        
        # Panes for some metric to track 1
        self.ax2 = self.fig.add_subplot(gs[0, -1])
        self.ax2.set_title('Metric 1')
        
        # Panes for some metric to track 2
        self.ax3 = self.fig.add_subplot(gs[1, -1])
        self.ax3.set_title('Metric 2')
        
        # Panes for some metric to track 3
        self.ax4 = self.fig.add_subplot(gs[2, -1])
        self.ax4.set_title('Metric 3')


river = River()
viz = Plotter(river, ships=None, sim_set=None)