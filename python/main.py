"""
This is the main function for the simulator. In here the important hyperparameters can
be set and the simulation is started

"""

# HYPERPARAMETERS

# Use Debugging
import matplotlib.animation as ani
import matplotlib.pyplot as plt
import visualization
from ships import Ships
from river import River
DEBUG = 0


# Initialize River
river = River()

# Number of vessels to spawn
N_SHIPS = 5

COORD_TYPE = "UTM"

# Length, width, and mass of each vessel.
# Y-Location of each ship during spawn
vessel_args = {
    "river": river,
    "num_ships": N_SHIPS,
    "ship_lengths": [100] * N_SHIPS,
    "ship_widths": [10] * N_SHIPS,
    "ship_mass": [100E3] * N_SHIPS,
    "y_location": [150] * N_SHIPS,
    "overtake_level": [0] * N_SHIPS,
    "direction": [1] * N_SHIPS,
    "spawn_dist": 300,
    "start_loc": 40_000}

# Init vessel class with pre-defined args
ships = Ships(**vessel_args)


# Granularity of timesteps
SIM_TIMESTEP = 2

# Which vessel is followed in the visualization
FOLLOWED_VESSEL = 0

if COORD_TYPE == "UTM":
    plotter = visualization.UTM_Plotter(
        river, ships, FOLLOWED_VESSEL, SIM_TIMESTEP)
else:
    plotter = visualization.Plotter(
        river, ships, FOLLOWED_VESSEL, SIM_TIMESTEP)


# Use Maplotlibs internal animation framework to speed up viz
def _update(i):

    water_depth = river.get_water_depth(ships)
    stream_vel = river.mean_stream_vel(ships)
    river_profile = river.get_river_profile(ships)

    for ship in range(N_SHIPS):
        ships.simulate_timestep(ship, river, SIM_TIMESTEP,
                                water_depth, river_profile, stream_vel)

    if COORD_TYPE == "UTM":

        for idx, ship in enumerate(plotter.vessel_obj):
            ship.set_xy(
                list(zip(*plotter.box_to_utm(ships.heading_box[idx]).exterior.xy)))

        x_utm, y_utm = river.get_utm_position(ships.x_location[plotter.followed_vessel],
                                              ships.y_location[plotter.followed_vessel])

        plotter.ax1.set_xlim(x_utm - 1000, x_utm + 1000)
        plotter.ax1.set_ylim(y_utm - 1000, y_utm + 1000)

    else:
        for idx, ship in enumerate(plotter.vessel_obj):
            ship.set_xy(list(zip(*ships.heading_box[idx].exterior.xy)))

        plotter.ax1.set_xlim(ships.x_location[plotter.followed_vessel] - 1000,
                             ships.x_location[plotter.followed_vessel]+1000)

    return plotter.vessel_obj

# This function uses a naive animation plot, however debugging is much easier
# while using it.


def _debug():
    while True:
        water_depth = river.get_water_depth(ships)
        stream_vel = river.mean_stream_vel(ships)
        river_profile = river.get_river_profile(ships)
        for ship in range(N_SHIPS):
            ships.simulate_timestep(
                ship, river, SIM_TIMESTEP, water_depth, river_profile, stream_vel)
        plotter.update()


if __name__ == "__main__":

    if DEBUG:
        _debug()
    else:
        simul = ani.FuncAnimation(
            plotter.fig, _update, blit=False, frames=200, interval=40)
        plt.show()
