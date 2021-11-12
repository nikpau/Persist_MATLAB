"""

Docstring again

"""
import matplotlib
matplotlib.use("TKAgg")
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.gridspec as gridspec
import matplotlib.animation as ani
import numpy as np
from shapely import geometry

COLOR = 'white'
matplotlib.rcParams['text.color'] = COLOR
matplotlib.rcParams['axes.labelcolor'] = COLOR
matplotlib.rcParams['xtick.color'] = COLOR
matplotlib.rcParams['ytick.color'] = COLOR

#from river import River
#from ships import Ships

class UTM_Plotter:

    def __init__(self, river, ships, followed_vessel: int, dT) -> None:
        
        print("Configuring UTM Plotter...")

        self.river = river
        self.ships = ships

       # Define colors:
        self.canvas_bg = "#1b4332"
        self.plot_bg = "#495057"

        self.followed_vessel = followed_vessel

        # Get number of spawned ships on the plane
        #self.num_ships = ships.num_ships

        self.fig = plt.figure()
        self.fig.patch.set_facecolor(self.canvas_bg)

        gs = self.fig.add_gridspec(3,3)

        # Pane for the map
        self.ax1 = self.fig.add_subplot(gs[:, :-1])
        self.ax1.set_title('Some Map')

        # Make a contour plot with river coords as x and y, and stream velocity as z (why +10?)
        self.ax1.contourf(river.point_coords[:,:,0], river.point_coords[:,:,1],river.stream_vel + 10,
        cmap = cm.winter_r)
        self.ax1.set_facecolor(self.plot_bg)

        # Plot Vessels
        # xy coords are unpacked and plotted simultaneously
        print("Placing vessels")
        self.vessel_obj = [plt.Polygon(list(zip(*self.box_to_utm(self.ships.heading_box[i]).exterior.xy)))\
             for i in range(self.ships.num_ships)]
        
        for patch in self.vessel_obj:
            patch.set_facecolor("black")
            self.ax1.add_patch(patch)

        # Set x and y limits to center around the selected vessel
        x_utm, y_utm = self.river.get_utm_position(self.ships.x_location[followed_vessel],
                                                    self.ships.y_location[followed_vessel])

        self.ax1.set_xlim(x_utm -500, x_utm + 500)
        self.ax1.set_ylim(y_utm -500, y_utm + 500)
        
        # Panes for some metric to track 1
        self.ax2 = self.fig.add_subplot(gs[0, -1])
        self.ax2.set_title('Metric 1')
        self.ax2.set_facecolor(self.plot_bg)
        
        # Panes for some metric to track 2
        self.ax3 = self.fig.add_subplot(gs[1, -1])
        self.ax3.set_title('Metric 2')
        self.ax3.set_facecolor(self.plot_bg)
        
        # Panes for some metric to track 3
        self.ax4 = self.fig.add_subplot(gs[2, -1])
        self.ax4.set_title('Metric 3')
        self.ax4.set_facecolor(self.plot_bg)

        self.fig.subplots_adjust(right = 0.98, left= 0.04, top= 0.96, bottom = 0.04)

    def update(self):

        self.ax1.clear()
        self.ax1.contourf(self.river.point_coords[:,:,0], self.river.point_coords[:,:,1],self.river.stream_vel + 10,
                            cmap = cm.winter_r)
        self.ax1.set_facecolor(self.plot_bg)
        for ship in range(self.ships.num_ships):
            self.ax1.fill(*self.box_to_utm(self.ships.heading_box[ship]).exterior.xy, "k")

        self.ax1.set_xlim(self.ships.x_location[self.followed_vessel] -1500,
                        self.ships.x_location[self.followed_vessel]+1500)
        # Set x and y limits to center around the selected vessel
        x_utm, y_utm = self.river.get_utm_position(self.ships.x_location[self.followed_vessel],
                                                    self.ships.y_location[self.followed_vessel])

        self.ax1.set_xlim(x_utm -500, x_utm + 500)
        self.ax1.set_ylim(y_utm -500, y_utm + 500)
        plt.pause(.001)


    def box_to_utm(self, box):

        x,y = box.exterior.xy

        new_x = np.zeros(4)
        new_y = np.zeros(4)

        for i in range(4):
            new_x[i], new_y[i] = self.river.get_utm_position(x[i],y[i])

        new_polygon = geometry.Polygon([*zip(new_x,new_y)])

        return new_polygon


class Plotter:
    def __init__(self, river, ships, followed_vessel: int, dT) -> None:
        

        print("Configuring Plotter...")

        self.river = river
        self.ships = ships

        self.followed_vessel = followed_vessel

        # Define colors:
        self.canvas_bg = "#1b4332"
        self.plot_bg = "#495057"

        # Get number of spawned ships on the plane
        #self.num_ships = ships.num_ships
        self.fig = plt.figure()
        self.fig.patch.set_facecolor(self.canvas_bg)

        gs = self.fig.add_gridspec(4,1)

        # Pane for the map
        self.ax1 = self.fig.add_subplot(gs[0, :])
        self.ax1.set_title('Some Map')

        # Make a contour plot with river coords as x and y, and stream velocity as z (why +10?)
        print("Initialize Mesh...")
        self.xx,self.yy = np.meshgrid(np.arange(0,520,20),
        np.array([n*20+1 for n in range(len(self.river.point_coords))]))
        self.ax1.contourf(self.yy,self.xx,river.stream_vel + 10, cmap = cm.winter_r)
        self.ax1.set_facecolor(self.plot_bg)

        # Plot Vessels
        # xy coords are unpacked and plotted simultaneously
        print("Placing vessels")
        self.vessel_obj = [plt.Polygon(list(zip(*self.ships.heading_box[i].exterior.xy)))\
             for i in range(self.ships.num_ships)]
        
        for patch in self.vessel_obj:
            patch.set_facecolor("black")
            self.ax1.add_patch(patch)

        # Set x limits to center around the selected vessel
        self.ax1.set_xlim(self.ships.x_location[followed_vessel] -2000,
                        self.ships.x_location[followed_vessel]+2000)
        
        # Panes for some metric to track 1
        self.ax2 = self.fig.add_subplot(gs[1, :])
        self.ax2.set_title('Metric 1')
        self.ax2.set_facecolor(self.plot_bg)
        
        # Panes for some metric to track 2
        self.ax3 = self.fig.add_subplot(gs[2, :])
        self.ax3.set_title('Metric 2')
        self.ax3.set_facecolor(self.plot_bg)
        
        # Panes for some metric to track 3
        self.ax4 = self.fig.add_subplot(gs[3,:])
        self.ax4.set_title('Metric 3')
        self.ax4.set_facecolor(self.plot_bg)

        # Adjust horizontal spacing between plots
        self.fig.subplots_adjust(hspace=0.25, right = 0.98, left=0.04, bottom = 0.04, top=0.96)