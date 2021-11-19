import gym
from gym import spaces
import numpy as np
import pandas as pd

"""
This is the main simulating environment:
   - There are n vessels placed on arbitary positions around the agent ship
   - The ships are spawned in the real river cooridnate system 
   - The agents gets rewarded for ...
"""

class RhineSim(gym.Env):
    def __init__(self, river_env, n_vessel = 5, v2v_dist = 300) -> None:

        # Init Gym Env
        super().__init__()

        # Number of vessels to be present at any time
        self.n_vessel = n_vessel

        # Vessel-to-vessel distance (Distance of between two non-agent vesssels)
        self.v2v_dist = v2v_dist

        # Make the rhine data globally available to the methods of the class
        self.river_env = river_env

        # Set the inital properties of the agent
        agent_properties = [
            "start_x_agent",
            "start_y_agent",
            "vx_agent",
            "vy_agent",
            "ax_agent"
            ,"ay_agent"
            ,"agent_lookahead" # Maximum foresight of agent. All events above this distance are irrelevant
            ]

        for property in agent_properties:
            setattr(self,property, 0.)

        # Set boundary conditions for the vessel behavior
        self.vx_max = 5
        self.vy_max = 5
        self.ax_max = 0
        self.ay_max = 0.01

        # Stepsize, current time step, max lenth of steps per episode
        self.dT = 1
        self.timestep = 0
        self.max_episode_steps = 2000

        



class JerkRiver:
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

# Use a box blur algorithm for the water_depth array to be smoothed
def box_blur(matrix: np.array) -> np.array:
    pass