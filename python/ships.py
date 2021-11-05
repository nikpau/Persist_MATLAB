"""
Arbitrary Docstring

"""
# from river import River

import numpy as np
from shapely import geometry, affinity

from river import River
from ships import Ships
import policies

class Ships:

    def __init__(self, 
                river, 
                num_ships, 
                ship_lengths, 
                ship_widths, 
                ship_mass,
                y_location,
                spawn_dist = 300,
                start_loc = 40000,
                ) -> None:
        
        assert isinstance(num_ships,int), "Only Integers allowed"

        #Set the number of ships to be spawned
        self.num_ships = num_ships

        self.ship_id = np.arange(self.num_ships)

        # Define the dimensions of the ships to be spawned
        assert len(ship_widths) == len(ship_lengths) == len(ship_mass) == num_ships, \
            "Dimensions must be specified for every ship"

        self.length = np.array(ship_lengths)
        self.width = np.array(ship_widths)
        self.eff_width = np.array(self.width) # ?? What is this?
        self.mass = np.array(ship_mass)

        # Absolute spawn point as river x-coordinate. (e.g. 30_000)
        self.start_location = start_loc

        # Distance between each spawned ship
        self.spawn_dist = spawn_dist

        # x positions for each ship
        self.x_location = self.start_location - np.linspace(0, self.num_ships * self.spawn_dist, self.num_ships)

        # Overtaking level per ship. (0 = no overtakling, 1 = overtaking, 2 = ??)
        self.overtaking_level = np.zeros(self.num_ships)

        # Direction each ship travels upon spawning (1 = upstream , -1 = downstream)
        self.direction = np.ones(self.num_ships)
        
        # Y location for each ship (must be an array with one entry per ship)
        self.y_location = np.array(y_location)

        # Velocity in x and y direction
        self.vx = 3 * self.direction
        self.vy = np.zeros(self.num_ships)

        # Acceleration in x and y direction
        self.ax = np.zeros(self.num_ships)
        self.ay = np.zeros(self.num_ships)

        # Set power per ship
        self.power = np.zeros(self.num_ships)

        # Squat
        self.squat = np.ones(self.num_ships)

        # Heading cf value 
        self.heading_cf = np.zeros(self.num_ships)

        # Heading angle in Radians
        self.heading_angle = np.zeros(self.num_ships)

        # Three points for circle calculation. These are num_ships' 2x3 matirces
        self.heading_radius_calc = np.zeros((self.num_ships,2,3))

        # List of objects holding the dimensions and alignment of each ship
        self.heading_box = [self.create_box(ID) for ID in range(self.num_ships)]

        # x and y coordinates (utm transformed)
        self.x_utm = np.zeros(self.num_ships)
        self.y_utm = np.zeros(self.num_ships) 

        # Vectors for the lateral and longitudinal control policy
        self.lat_con_pol = []
        self.long_con_pol = []

        # TODO: Init control policies for each ship (Ships.m:105-108)

    # Create a Polygon for each ship
    def create_box(self, ID):

        l = self.length[ID]
        w = self.width[ID]

        # Create Polygon with edges of the ship
        polygon = geometry.Polygon(((-w/2,l/2),(w/2,l/2),(w/2,-l/2),(-w/2,-l/2)))

        # Calculate the heading angle of the ship
        self.heading_angle[ID] = self.heading_angle[ID] + np.arctan(self.vy[ID]/self.vx[ID]) / (2 * np.pi * 360)

        # Rotate the polygon around the heading angle
        rotated_polygon = affinity.rotate(polygon, self.heading_angle[ID])

        # Move the unity polygon to the respective point on the river
        box = affinity.translate(rotated_polygon, self.y_location[ID], self.x_location[ID])

        return box

    # Find the ship ahead or behind with the minimum distance to the agent
    def find_ships(self, ID, num, overtake_level, mode = "ahead"):

        if mode == "ahead":
            dist_with_lengths = self.direction[ID] * (self.x_location + self.direction[ID]*self.length / 2 - \
                (self.x_location[ID] - self.direction[ID] * self.length[ID] / 2))
            dist = self.direction[ID] * (self.x_location - self.x_location[ID])

        elif mode in ["behind","oncoming"]:
            dist_with_lengths = self.direction[ID] * (self.x_location[ID] + self.direction[ID]*self.length[ID] / 2 - \
                (self.x_location - self.direction[ID] * self.length[ID] / 2))
            dist = self.direction[ID] * (self.x_location[ID] - self.x_location)

        # Set distances to other ships to 'Inf' if they are the agent, have different overtaking levels,
        # or drive in different directions
        if mode == "oncoming":
            subsetter = np.where((dist_with_lengths <= 0) | \
                            (self.direction == self.direction[ID]) |\
                            (self.ship_id == ID) | \
                            (self.overtaking_level != overtake_level))
        else:
            subsetter = np.where((dist_with_lengths <= 0) | \
                            (self.direction != self.direction[ID]) |\
                            (self.ship_id == ID) | \
                            (self.overtaking_level != overtake_level))

        dist[subsetter] = float("Inf")

        # Sort the distances in ascending order
        sorted_dist = np.sort(dist)

        # Get the first num minimum distances
        min_dist = sorted_dist[1:(num + 1)]

        indices = np.full(len(min_dist), float("Inf"))

        # Set the indices to where the min dist appears in the original dist array
        for i in range(len(min_dist)):
            if min_dist[i] == float("Inf"):
                indices[i] = float("Inf")
            else:
                indices[i] = np.where(dist == min_dist[i])[0][0]

        return (indices, min_dist)

    # TODO: Implement Ships.m:232
    # We compute some heading based on a cf value?
    def compute_heading_from_cf(self, ID):

        current_points = self.heading_radius_calc[ID]
        added_coords = np.c_[current_points, [self.x_utm,self.y_utm]]

        if current_points.all() == 0:
            heading_angle = 0.

        else:
            # Remove the first of the three points in order to achieve a rolling window
            added_coords = np.delete(added_coords,0,1)
            r, xy = self.fit_circle(added_coords)
            heading_angle = np.arcsin(self.length[ID] * (self.heading_cf[ID] - 0.5) / r)

        # Calculate whether the center_point is
        # on the right or left of the line spanned by two points of
        # our added_coords list. This is most easily achieved using a cross product
        if np.cross(added_coords[1] - added_coords[0], xy - added_coords[0]) > 0:
            # Point is on the left
            pass



    # Find the unique circle inscribed by three points using the inscribed angle theorem
    @staticmethod
    def fit_circle(points: np.array):

        x = points[0,:]
        y = points[1,:]

        x_center = (x[0]**2 - x[1]**2 + y[0]**2 - y[1]**2) / 2*(x[0] - x[1])
        y_center = (y[0]**2 - y[2]**2 + x[0]**2 - x[2]**2) / 2*(y[0] - y[2])

        radius = np.sqrt((x[0] - x_center)**2 + (y[0] - y_center)**2)

        return radius, (x_center,y_center)


river = River()
ships = Ships(river=river, num_ships=5, ship_lengths= [100]*5, ship_widths= [10]*5,ship_mass=[1e5]*5,y_location=[150]*5)
lat = policies.LatConPol(30,1)
f = lat.compute_obs(ships, river,2)