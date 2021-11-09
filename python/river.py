"""
This will be a Docstring in the future

"""

from numpy.core.numeric import indices
import pandas as pd
import numpy as np 


class River:
    def __init__(self) -> None:
        
        # Number of evenly spaced gridpoints for which stream velocity and water depth data is available. 
        # All Gridpoints are BASEPOINT_DIST meters apart to span a maximum rhine width of 500 meters 
        self.GRIDPOINTS_PER_WIDTH = 26
        self.BASEPOINT_DIST = 20

        self.min_water_under_keel = 0.2

        # Coordninates for the river loctation
        self.point_coords = pd.read_csv("python/data/MW_grid_500m_655_852_UTM32_1.txt", sep=" ", skiprows=1)
        self.point_coords = self.point_coords.iloc[:,1:3] # Only longitudinal and lateral datapoints

        # Data for water depth and stream velocity
        self.point_dataset = pd.read_csv("python/data/MW_grid_500m_655_852_UTM32_2.txt", sep=" ")
        self.point_dataset = self.point_dataset.iloc[:,0:6]

        # Init vectors for water depth and stream velocity
        self.water_depth = self.stream_vel = np.zeros((len(self.point_coords)))

        self.stream_vel = np.sqrt(np.add(self.point_dataset.iloc[:,0].values**2, self.point_dataset.iloc[:,1].values**2))
        self.water_depth = self.point_dataset.iloc[:,2].values

        # Reshape stream velocity and water depth to have GRIDPOINTS_PER_WIDTH columns
        self.stream_vel = self.stream_vel.reshape((-1, self.GRIDPOINTS_PER_WIDTH))
        self.water_depth = self.water_depth.reshape((-1, self.GRIDPOINTS_PER_WIDTH))

        self.point_coords = self.point_coords.values.reshape((-1,self.GRIDPOINTS_PER_WIDTH,2))        
        self.point_dataset = self.point_dataset.values.reshape((-1,self.GRIDPOINTS_PER_WIDTH,6))

        # Create upper and lower basepoints

        self.length = self.BASEPOINT_DIST * len(self.point_coords)

        mean_base_point_dist =  500

        # Create equally spaces points for the entire river length, and add some random noise to them.
        self.base_points_l = np.arange(start=0,stop=self.length,step=mean_base_point_dist)
        self.base_points_l = np.add(self.base_points_l,np.random.normal(loc=0,scale=10,size=len(self.base_points_l)))

        self.base_points_u = np.arange(start=0,stop=self.length,step=mean_base_point_dist)
        self.base_points_u = np.add(self.base_points_u, np.random.normal(loc=0,scale=10,size=len(self.base_points_u)))

        # Add second dimension with all zeros
        zeros = np.zeros(len(self.base_points_u))
        self.base_points_l = np.append(self.base_points_l, zeros).reshape(2,-1)
        self.base_points_u = np.append(self.base_points_u, zeros).reshape(2,-1)

    def get_utm_position(self, x, y):

        x_log = int(x/self.BASEPOINT_DIST) # Casting to int, i.e. flooring
        x_alpha = x/self.BASEPOINT_DIST - x_log
        
        y_log = int(y/self.BASEPOINT_DIST) + 1 # Casting to int, i.e. flooring
        y_alpha = y/self.BASEPOINT_DIST + 1 - y_log

        # Base coordinates after utm transform
        base_utm = self.point_coords[x_log,y_log,:]

        # Difference between points in x direction
        diff_x = self.point_coords[x_log+1,y_log,:] - base_utm

        base_utm += x_alpha * diff_x

        # Difference between points in y direction
        diff_y = self.point_coords[x_log,y_log+1,:] - base_utm

        x_utm = base_utm[0] + y_alpha * diff_y[0]
        y_utm = base_utm[1] + y_alpha * diff_y[1]

        return [x_utm, y_utm]

        
    # This is the same function as in River.m:148-157 but vectorized. 
    # Therefore mean stream velocities for every vessel are calculated at once.
    # This means the function only gets called once per timestep and not once per ship per timestep
    def mean_stream_vel(self, ships):

        min_y_log = np.round((ships.y_location - ships.eff_width / 2) / self.BASEPOINT_DIST).astype(int)
        max_y_log = np.round((ships.y_location + ships.eff_width / 2) / self.BASEPOINT_DIST).astype(int)

        min_x_log = np.round((ships.x_location - ships.eff_width / 2) / self.BASEPOINT_DIST).astype(int) + 10
        max_x_log = np.round((ships.x_location + ships.eff_width / 2) / self.BASEPOINT_DIST).astype(int) + 10

        # Define empty array for all mean stream velocities to be written into
        v_stream = np.zeros(len(min_y_log))
        for i in range(len(v_stream)):
            v_stream[i] = ships.direction[i] *\
                 np.mean(self.stream_vel[min_x_log[i]:max_x_log[i]+1, min_y_log[i]:max_y_log[i]+1])

        return v_stream

    # Vectorized function for receiving the current water depth per vessel (orig: River.m:169)
    def get_water_depth(self,ships):

        min_x_log = np.round((ships.x_location - ships.length / 2) / self.BASEPOINT_DIST).astype(int) - 10
        max_x_log = np.round((ships.x_location + ships.length / 2) / self.BASEPOINT_DIST).astype(int) + 10

        min_y_log = np.floor(ships.y_location / self.BASEPOINT_DIST).astype(int)
        max_y_log = np.ceil(ships.y_location / self.BASEPOINT_DIST).astype(int)

        depth = np.zeros(ships.num_ships)

        for i in range(len(depth)):
            depth[i] = np.mean(np.mean(self.water_depth[min_x_log[i]:max_x_log[i]+1,min_y_log[i]:max_y_log[i]+1]))

        return depth

    # Find the river profile for every ship (orig River.m:159)
    def get_river_profile(self,ships):

        min_x_log = np.round((ships.x_location - ships.length / 2) / self.BASEPOINT_DIST).astype(int) - 10
        max_x_log = np.round((ships.x_location + ships.length / 2) / self.BASEPOINT_DIST).astype(int) + 10

        river_profile = np.zeros(ships.num_ships)

        for i in range(len(river_profile)):
            river_profile[i] = np.mean(np.trapz(self.water_depth[min_x_log[i]:max_x_log[i]+1,:])) * self.BASEPOINT_DIST

        return river_profile

    # Vectorized version of the free space under keel function in River.m:178
    def detect_free_space(self, ID, ships, look_ahead_dist, num_points):

        look_ahead_dist = look_ahead_dist / self.BASEPOINT_DIST

        # Get x position of vessel by subtracting its length from its current x coordinate
        pos_x = np.round((ships.x_location[ID] - ships.length[ID] / 2) / self.BASEPOINT_DIST)

        # Get vessel squat values
        ssq = ships.squat[ID]

        # Split the lookahead into num_points blocks for lateral control policy
        length_block = look_ahead_dist / num_points


        if ships.direction[ID] == 1:
            start = np.round((np.floor(pos_x/length_block) - 1) * length_block)
        else:
            start = np.round((np.ceil(pos_x/length_block) + 1) * length_block) - look_ahead_dist + 1
        
        # This returns a list with num_ships arrays containing the respective water depths [water_depth1, water_depth2,...]
        water_depth_array  = self.water_depth[int(start)-1:(int(start) + int(look_ahead_dist) -1), :]

        # List containing logical arrays per ship with 1 being coast and 0 being free space
        free_space = water_depth_array - (ssq + self.min_water_under_keel) < 0

        # Extent the coast indicator points (ones) per length-block over the length of the vessel
        # to avoid it to reeve into the coast because there is clearance in front
        #
        # This is achieved by duplicating the freespace matrix and adding one column of zeros 
        # to the beginning of the first and to the end of the second. This moves the first matrix one step in x direction
        # Now both matrices are iteratively XOR-ed to carry over high points (ones) from one step to the next, 
        # thereby extending high points of the matix over time.
        for i in range(int(np.round(ships.length[ID] / self.BASEPOINT_DIST))):
            a = np.r_[free_space, [np.zeros(free_space.shape[1])]]
            b = np.r_[[np.zeros(free_space.shape[1])], free_space]
            free_space = np.array(a, dtype=bool) | np.array(b, dtype=bool)

        # Find the first points in the logic freespace array that belong to the coast
        # Freespace array: 26xN
        # |1 1 1 1 1 1 1 1 1 1| > y1
        # |0 1 0 0 1 0 0 0 0 1| > y ..
        # |0 0 0 0 0 0 0 1 0 0| > y ..
        # |1 1 1 1 1 1 1 1 1 1| > y_26
        # - dir of vessel (N) ->
        # 
        # Here the min and max y-indices that contain a transision from 0 to one are returned
        # to determine where the vessel can safely drive
        #
        lower_poly = np.zeros(free_space.shape[0])
        upper_poly = np.zeros(free_space.shape[0])

        for i in range(free_space.shape[0]):
            added_ones = np.insert(free_space[i], [0,len(free_space[i])], 1).astype(np.float32)
            diff_of_vals = np.abs(np.diff(added_ones).astype(np.float32))
            non_zero_indices = np.where((diff_of_vals != 0))
            diff_of_indices = np.diff(non_zero_indices)
            max_index = np.argmax(diff_of_indices)


            lower_poly[i] = non_zero_indices[0][max_index]
            upper_poly[i] = non_zero_indices[0][max_index + 1]

        # From here on the results of the matlab implementation differ from this one due to zero indexing.
        # The linspacer starts at one as it does in the original file River.m:218.
        # However, the first one-index is used to index the first element of the upper_poly vector in l. 203.
        # In matlab this works due to one-indexing. In python it does not. Therefore the first element in 
        # the upper_poly array is omitted. If the linspacer would start at zero, the breakpoints of the spacer
        # would differ from the original. Difference is minor, therefore no further adjustments are made.
        blocks_upper = np.round(np.linspace(1,len(upper_poly),num_points + 1)).astype(int)
        blocks_lower = np.round(np.linspace(1,len(upper_poly),num_points + 1)).astype(int)

        upper_points = np.zeros((num_points,2))
        lower_points = np.zeros((num_points,2))

        for i in range(num_points):
            
            current_upper_range = upper_poly[blocks_upper[i]:blocks_upper[i+1]]
            current_lower_range = lower_poly[blocks_lower[i]:blocks_lower[i+1]]

            upper_min = np.min(current_upper_range)
            upper_ind = np.argmin(current_upper_range) 
            upper_points[i,:] = upper_ind + blocks_upper[i], upper_min

            lower_max = np.max(current_lower_range)
            lower_ind = np.argmax(current_lower_range)
            lower_points[i,:] = lower_ind + blocks_lower[i], lower_max

        upper_points[:,0] = (upper_points[:,0] + start -1) * self.BASEPOINT_DIST
        upper_points[:,1] = (upper_points[:,1]) * self.BASEPOINT_DIST
        lower_points[:,0] = (lower_points[:,0] + start -1) * self.BASEPOINT_DIST
        lower_points[:,1] = (lower_points[:,1]) * self.BASEPOINT_DIST

        return (lower_points, upper_points)
