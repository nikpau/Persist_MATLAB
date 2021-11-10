"""
This is the main function for the simulator. In here the important hyperparameters can
be set and the simulation is started

"""

# HYPERPARAMETERS

from river import River
from ships import Ships


N_SHIPS = 5

# Length, width, and mass of each vessel. 
SHIP_LENGTH = [100] * N_SHIPS
SHIP_WIDTH = [20] * N_SHIPS
SHIP_MASS = [100E3] * N_SHIPS

# Y-Location of each ship during spawn
Y_LOCAION = [150] * N_SHIPS

# Granularity of timesteps
SIM_TIMESTEP = 1

# Which vessel is followed in the visualization
FOLLOW_VESSEL = 0

# Initialize River and ships
river = River()
ships = Ships(river, N_SHIPS,SHIP_LENGTH,SHIP_WIDTH,SHIP_MASS,Y_LOCAION)

def main():
    # Init simulation timesteps
    timestep = 0

    while True:

        timestep += 1

        water_depth = river.get_water_depth(ships)
        stream_vel = river.mean_stream_vel(ships)
        river_profile = river.get_river_profile(ships)

        for ship in range(N_SHIPS):
            ships.simulate_timestep(ship, river,SIM_TIMESTEP,water_depth,river_profile,stream_vel)
