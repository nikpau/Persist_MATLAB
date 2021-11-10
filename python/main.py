"""
This is the main function for the simulator. In here the important hyperparameters can
be set and the simulation is started

"""

# HYPERPARAMETERS

from river import River
from ships import Ships
import visualization

# Initialize River 
river = River()

# Number of vessels to spawn
N_SHIPS = 5

# Length, width, and mass of each vessel. 
# Y-Location of each ship during spawn
vessel_args = {
        "river":river,
        "num_ships": N_SHIPS,
        "ship_lengths":[100] * N_SHIPS,
        "ship_widths":[10] * N_SHIPS,
        "ship_mass":[100E3] * N_SHIPS,
        "y_location": [150] * N_SHIPS,
        "overtake_level": [0] * N_SHIPS,
        "direction": [1] * N_SHIPS,
        "spawn_dist": 300,
        "start_loc": 40_000}

# Init vessel class with pre-defined args
ships = Ships(**vessel_args)


# Granularity of timesteps
SIM_TIMESTEP = 1

# Which vessel is followed in the visualization
FOLLOWED_VESSEL = 0


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
