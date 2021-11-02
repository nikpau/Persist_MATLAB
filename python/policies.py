
"""
In here, the lateral and longitudinal contol policies are definied

"""

from ships import Ships
from river import River

import numpy as np

class LatConPol:

    def __init__(self, ID) -> None:
        
        self.id = ID

        self.upper_speed_bound = 4
        self.upper_ttc_bound = 1000
        self.lower_ttc_bound = -1000
        self.upper_river_bound = 500
        self.upper_acc_bound = 0.01
        self.safety_dist = 0
        self.num_features = 14
        self.look_ahead_dist = 4500

        # Number of breakpoints the lookahead distance gets divided into
        self.n_breakpoints = 7

    def compute_obs(self, ships, river, ships_to_consider):

        ships_to_consider += 1 # Python has open intervals

        # Init time to collision matricies
        self.TTC_ship_ahead = np.zeros((ships_to_consider,2))
        self.TTC_y_ship_ahead = np.zeros((ships_to_consider,2))
        
        self.TTC_ship_behind = np.zeros((ships_to_consider,2))
        self.TTC_y_ship_behind = np.zeros((ships_to_consider,2))

        self.TTC_ship_oncoming = np.zeros((ships_to_consider,2))
        self.TTC_y_ship_oncoming = np.zeros((ships_to_consider,2))

        # Extract ship ID and direction
        ID = self.id
        direction = ships.direction[ID]


        (self.base_points_lower,self.base_points_upper) = river.detect_free_space(ID, ships, self.look_ahead_dist, self.n_breakpoints)

        # Calculate time to collision vector for the upper and lower coast line
        self.TTC_base_l = direction * (self.base_points_lower[:,0] - ships.x_location[ID]) / np.abs(ships.vx[ID])
        self.TTC_base_l = np.tile(self.TTC_base_l, (2,1))
        self.y_base_l = self.base_points_lower[:,1]

        self.TTC_base_u = direction * (self.base_points_upper[:,0] - ships.x_location[ID]) / np.abs(ships.vx[ID])
        self.y_base_u = self.base_points_upper[:,1]

        # Ship ahead observation
        overtake_range_ahead = np.arange(start=0 ,stop=ships.overtaking_level[ID])
        overtake_range_behind = np.arange(start=ships.overtaking_level[ID]+1, stop=4)
        overtake_range_oncoming = np.arange(start=0 ,stop=4)

        (ahead_ids, _) = ships.find_ships(ID, ships_to_consider, overtake_range_ahead, mode = "ahead")
        (behind_ids, _) = ships.find_ships(ID, ships_to_consider, overtake_range_behind, mode = "behind")
        (oncoming_ids, _) = ships.find_ships(ID, ships_to_consider, overtake_range_oncoming, mode = "oncoming")
        
        # Exterior coordinates of the agent ship
        ext_coords_agent = np.flip(np.array(ships.heading_box[ID].exterior.coords),1)

        for i in range(ships_to_consider):

            # ------------------------
            # SHIPS AHEAD
            # ------------------------

            # TODO Double check this. There are some errors in here for sure


            # Ship ahead is either not considered or faster than the agent,
            # therefore time to collision is not needed, i.e. infinite.
            if ahead_ids[i] == float("Inf"):
                self.TTC_ship_ahead[i,:] = float("Inf"), float("Inf")
                self.TTC_y_ship_ahead[i,:] = 0., 0.
            
            elif (np.abs(ships.vx[ID]) < np.abs(ships.vx[ahead_ids[i]])):

                self.TTC_ship_ahead[i,:] = float("Inf"), float("Inf")
                self.TTC_y_ship_ahead[i,:] = 0., 0.

            # Ship ahead slower and considered
            else:
                ext_coords_ship_ahead = np.flip(np.array(ships.heading_box[ahead_ids[i]].exterior.coords),1)
                if direction == 1:

                    # X - Difference
                    x_difference = [
                        ext_coords_ship_ahead[3,0] - ext_coords_agent[1,0],
                        ext_coords_ship_ahead[0,0] - ext_coords_agent[2,0]
                    ]
                    

                    # Y - Difference
                    y_difference = [
                        ext_coords_ship_ahead[3,1] - ext_coords_agent[1,1],
                        ext_coords_ship_ahead[0,1] - ext_coords_agent[2,1]
                    ]

                
                else:

                    # X - Difference
                    x_difference = [
                        ext_coords_agent[3,0] - ext_coords_ship_ahead[1,0],
                        ext_coords_agent[2,0] - ext_coords_ship_ahead[2,0]
                    ]
                    

                    # Y - Difference
                    y_difference = [
                        ext_coords_agent[3,1] - ext_coords_ship_ahead[1,1],
                        ext_coords_agent[2,1] - ext_coords_ship_ahead[0,1]
                    ]

                # Calc x difference of points divided by difference of speeds 
                # to get time to collision
                self.TTC_ship_ahead[i,:] = x_difference / (np.abs(ships.vx[ID] - ships.vx[ahead_ids[i]]))

                self.TTC_y_ship_ahead[i,:] = y_difference

            # ----------------------
            # SHIPS BEHIND
            # ----------------------

            # Ship behind is either not considered or slower than the agent,
            # therefore time to collision is not needed, i.e. infinite.
            if behind_ids[i] == float("Inf"):
                self.TTC_ship_behind[i,:] = float("Inf"), float("Inf")
                self.TTC_y_ship_behind[i,:] = 0., 0.

            elif (np.abs(ships.vx[ID]) > np.abs(ships.vx[ahead_ids[i]])):

                self.TTC_ship_behind[i,:] = float("Inf"), float("Inf")
                self.TTC_y_ship_behind[i,:] = 0., 0.

            else:
                ext_coords_ship_behind = np.flip(np.array(ships.heading_box[behind_ids[i]].exterior.coords),1)
                if direction == 1:

                    # X - Difference
                    x_difference = [
                        ext_coords_agent[3,0] - ext_coords_ship_behind[1,0],
                        ext_coords_agent[0,0] - ext_coords_ship_behind[2,0]
                    ]
                    

                    # Y - Difference
                    y_difference = [
                        ext_coords_agent[3,1] - ext_coords_ship_behind[1,1],
                        ext_coords_agent[0,1] - ext_coords_ship_behind[2,1]
                    ]

                
                else:

                    # X - Difference
                    x_difference = [
                        ext_coords_ship_behind[3,0] - ext_coords_agent[1,0],
                        ext_coords_ship_behind[0,0] - ext_coords_agent[2,0]
                    ]
                    

                    # Y - Difference
                    y_difference = [
                        ext_coords_ship_behind[3,1] - ext_coords_agent[1,1],
                        ext_coords_ship_behind[0,1] - ext_coords_agent[2,1]
                    ]

                # Calc x difference of points divided by difference of speeds 
                # to get time to collision
                self.TTC_ship_behind[i,:] = x_difference / (np.abs(ships.vx[ID] - ships.vx[behind_ids[i]]))

                self.TTC_y_ship_behind[i,:] = y_difference

            # ----------------------
            # SHIPS ONCOMING
            # ----------------------

            # Ship oncoming is not considered,
            # therefore time to collision is not needed, i.e. infinite.
            if oncoming_ids[i] == float("Inf"):

                self.TTC_ship_oncoming[i,:] = float("Inf"), float("Inf")
                self.TTC_y_ship_oncoming[i,:] = 0., 0.

            else:
                ext_coords_ship_oncoming = np.flip(np.array(ships.heading_box[oncoming_ids[i]].exterior.coords),1)
                if direction == 1:

                    # X - Difference
                    x_difference = [
                        ext_coords_ship_oncoming[2,0] - ext_coords_agent[0,0],
                        ext_coords_ship_oncoming[1,0] - ext_coords_agent[3,0]
                    ]
                    

                    # Y - Difference
                    y_difference = [
                        ext_coords_ship_oncoming[2,1] - ext_coords_agent[0,1],
                        ext_coords_ship_oncoming[1,1] - ext_coords_agent[3,1]
                    ]

                
                else:

                    # X - Difference
                    x_difference = [
                        ext_coords_agent[2,0] - ext_coords_ship_oncoming[1,0],
                        ext_coords_agent[1,0] - ext_coords_ship_oncoming[2,0]
                    ]
                    

                    # Y - Difference
                    y_difference = [
                        ext_coords_agent[2,1] - ext_coords_ship_oncoming[1,1],
                        ext_coords_agent[1,1] - ext_coords_ship_oncoming[2,1]
                    ]

                # Calc x difference of points divided by difference of speeds 
                # to get time to collision
                self.TTC_ship_behind[i,:] = x_difference / (np.abs(ships.vx[ID] - ships.vx[oncoming_ids[i]]))

                self.TTC_y_ship_behind[i,:] = y_difference

        if direction == 1:
            self.TTC_vector_up, self.POS_vector_up = self.finalize_ttc_vector(
                dir="upper",
                n_entries=self.num_features/2,
                coast_ttc=self.TTC_base_u, 
                vessel_ttc=np.concatenate([self.TTC_ship_behind,self.TTC_ship_oncoming]),
                coast_pos=self.y_base_u,
                vessel_pos=np.concatenate([self.TTC_y_ship_behind,self.TTC_y_ship_oncoming])
            )
            
            self.TTC_vector_low, self.POS_vector_low = self.finalize_ttc_vector(
                dir = "lower",
                n_entries= self.num_features/2,
                coast_ttc=self.TTC_base_l,
                vessel_ttc=self.TTC_ship_ahead,
                coast_pos=self.y_base_l,
                vessel_pos=self.TTC_y_ship_ahead
            )
        else:
            self.TTC_vector_up, self.POS_vector_up = self.finalize_ttc_vector(
                dir="upper",
                n_entries=self.num_features/2,
                coast_ttc=self.TTC_base_u, 
                vessel_ttc=self.TTC_ship_ahead,
                coast_pos=self.y_base_u,
                vessel_pos=self.TTC_y_ship_ahead
            )
            
            self.TTC_vector_low, self.POS_vector_low = self.finalize_ttc_vector(
                dir = "lower",
                n_entries= self.num_features/2,
                coast_ttc=self.TTC_base_l,
                vessel_ttc=np.concatenate([self.TTC_ship_behind,self.TTC_ship_oncoming]),
                coast_pos=self.y_base_l,
                vessel_pos=np.concatenate([self.TTC_y_ship_behind,self.TTC_y_ship_oncoming])
            )

        self.obs = np.array([
            ships.ay[ID] / self.upper_acc_bound,
            ships.vy[ID] / self.upper_speed_bound,
            self.TTC_vector_up / self.upper_ttc_bound,
            self.TTC_vector_low / self.upper_ttc_bound,
            np.tanh(self.POS_vector_up / 300),
            np.tanh(self.POS_vector_low / 300)
        ])

        return self.obs


    # In here, the vessels and coastline are being melted together to form a homogenous
    # coast line. This way the NN only gets a single TTC vector as input. 
    # As we melt the other vessels into the coastline we do not need the NN to decide between
    # vessel or coast.
    def finalize_ttc_vector(self, dir, n_entries, coast_ttc, vessel_ttc, coast_pos, vessel_pos):

        assert dir in ["upper","lower"],"Can only calulate upper or lower coastline"

        pos = np.concatenate((coast_pos, vessel_pos))
        ttc = np.concatenate((coast_ttc,vessel_ttc))

        # Sort the ttc and position vector in ascending or descending order depending on
        # the direction of the vessel
        if dir == "upper":
            pos = pos[pos[:,0].argsort()] # Sort in ascending order based on the first col
            ttc = ttc[ttc[:,0].argsort()]
        else:
            pos = pos[pos[:,0].argsort()[::-1]] # Sort in descending order based on the first col
            ttc = ttc[ttc[:,0].argsort()[::-1]]

        to_del = np.zeros((pos.shape[0],2))

        for i in range(pos.shape[0]-1):
            ttc_low = ttc[i,0]
            ttc_high = ttc[i,1]

            # Hidden tuple
            if ttc[i+1,0] > ttc_low and ttc[i+1,1] < ttc_high:

                ttc[i+1,0:2] = ttc_low, ttc_high
                to_del[i+1,0:2] = 1,1
            
            # Lower tuple
            elif ttc[i+1,0] < ttc_low and ttc[i+1,1] < ttc_high and ttc[i+1,1] > ttc_low:

                ttc[i+1,1] = ttc_low
                to_del[i+1,1] = 1

            # Upper tuple
            elif ttc[i+1,0] > ttc_low and ttc[i+1,1] > ttc_high and ttc[i+1,0] < ttc_low:

                ttc[i+1,0] = ttc_high
                to_del[i+1,0] = 1
            
        ttc = np.concatenate([ttc[:,0], ttc[:,1]])
        pos = np.concatenate([pos[:,0], pos[:,1]])
        to_del = np.concatenate([to_del[:,0], to_del[:,1]])

        # Delete all entries form the ttc and pos vectors according to the deletion vector
        ttc = [value for (value, bool) in zip(ttc,to_del) if bool]
        pos = [value for (value, bool) in zip(pos,to_del) if bool]
            
        # Sort both vectors based on the ttc vector in ascending order
        sort_indices = np.argsort(ttc)
        ttc = ttc[sort_indices]
        pos = pos[sort_indices]

        # Delete all double entries
        ttc = np.unique(ttc)
        pos = np.unique(pos)

        while ttc[1] < 0:
            pos = np.delete(pos,0)
            ttc = np.delete(ttc,0)

        if len(ttc) < n_entries:
            while len(ttc) < n_entries:
                np.append(ttc, self.upper_ttc_bound)
                np.append(pos, 0)
        
        pos = pos - self.safety_dist if dir == "upper" else pos + self.safety_dist

        TTC = np.clip(ttc,self.lower_ttc_bound,self.upper_ttc_bound)
        POS = np.clip(pos,-self.upper_river_bound,self.upper_river_bound)

        return [TTC,POS]


# Longitudinal control policy
class LonConPol:
    def __init__(self, ID, ships, river) -> None:

        self.id = ID

        # Make input args available in whole class
        self.river = river
        self.ships = ships

        # Init constants
        self.max_power = 1e6
        self.upper_speed_bound = 8
        self.max_dist_to_consider = 300
        self.v_str_scale = 4 # TODO: Whats this?
        self.a_mean = 1500
        self.a_scale = 800
        self.h_mean = 5.5
        self.h_scale = 3


        # Receive the current stream velocities for each vessel
        self.obs_stream_vel = river.mean_stream_vel(ships)

        # TODO What is this?
        self.wid = 0

    # Original in LonConPol.m:32
    # Here, the water_depth, stream velocity and river profile come from an external vectorized function
    # Therefore, we do not need to calculate the values for each vessel but can just read them which is O(1).
    def compute_obs(self, water_depth, stream_vel, river_profile):

        # Get the ID and the distance of the closest ship ahead
        ahead_ID, self.dist_ship_ahead = self.ships.find_ship(self.id, 1, self.ships.overtake_leve[self.id], mode = "ahead")

        # Get vertical velocity for the ship with minimal distance ahead.
        # If there is none in range, set the speed to the maximum possible value
        if ahead_ID is not float("Inf"):
            v_ahead = self.ships.vx[ahead_ID]
            self.dist_ship_ahead = self.dist_ship_ahead - self.ships.length[self.id]/2 - self.ships.length[ahead_ID]/2
        
        else:
            v_ahead = self.upper_speed_bound
            self.dist_ship_ahead = self.max_dist_to_consider

        # Clip the maximum distance to the maximum considerable distance
        self.dist_ship_ahead = np.clip(self.dist_ship_ahead,0,self.max_dist_to_consider) 

        # Import river properties
        self.obs_stream_vel = stream_vel
        self.river_profile = river_profile
        self.water_depth = water_depth

        self.obs = np.array(
            [
                self.ships.direction[self.id] * self.ships.vx[self.id] / self.upper_speed_bound,
                self.dist_ship_ahead / self.max_dist_to_consider,
                self.ships.direction[self.id] * (v_ahead - self.ships.vx[self.id]) / self.upper_speed_bound,
                self.ships.power[self.id] / self.max_power,
                self.obs_stream_vel / self.v_str_scale,
                (self.river_profile - self.a_mean) / self.a_scale,
                (self.water_depth - self.h_mean) / self.h_scale
            ]
        )

    # TODO What is ACC? Acceleration? Anyway, it is calculated here:
    def compute_acc(self, water_depth, stream_vel, river_profile, TIMESTEP):

        prev_stream_vel = self.obs_stream_vel

        # Import river properties
        self.obs_stream_vel = stream_vel
        self.river_profile = river_profile
        self.water_depth = water_depth

        # Set dynamic variables
        h = self.water_depth
        A = self.river_profile

        # absolute x-velocity of vessel
        v = np.abs(self.ships.vx[self.id])
        
        # Velocity relative to the river
        v_rel = v - self.obs_stream_vel

        # Constants for calculation
        RHO = 1000 # Water density
        GAMMA = 0.667 # Block coefficent
        LS = self.ships.length[self.id]
        BS = self.ships.width[self.id]
        TS = 2.8 # Draught + Additional Draught due to speed ("Abladetiefe")
        TSMAX = TS # ???
        AS = BS * TS
        TSAVG = 0.85 # Average draught

        MS = 695285 # TODO What is this

        CK = 0.8 # TODO What is this?
        KS = 0.2 # Rauigkeit Sohle [m]

        KSS = 0.0003 # Rauigkeit Unterschiff
        CSTR = 0.78 # What is this?
        ALPHA = 1.5E-4 # Also this?

        N = A / AS # And this
        VR = -v/(N - 1) # And this?
        MH = self.hydrodyn_mass(A, AS, MS)

        # Calculate the vessel squat for the next time step (after Barras)
        new_squat = 0.217 * N**(-0.76) * GAMMA * v_rel**2

        # Vessel resistane (Schiffswiderstand)
        V_REFF = -v_rel*(1 / (A / (BS * (TS + new_squat)) - 1))
        CW = 0.1 + 0.2 * (TSMAX / h)**2

        W1 = 0.5 * CW * CK * RHO * AS * (v_rel -V_REFF)**2

        lambda_s = (1.89 + 1.62 * np.log10(h/KSS))**(-2.5)

        W2 = 0.5 * CW * CK * RHO * lambda_s * BS * LS * (1 + 2 * TS / BS) * GAMMA * (v_rel - V_REFF)

        IDIR = -np.sign(self.obs_stream_vel - V_REFF)
        W3 = 0.023 * IDIR * CSTR * RHO * (KS * new_squat / (4 * h))**(1/3) * LS * BS * (TS / h) \
            * (1 + 2 * TS / BS) * GAMMA * (self.obs_stream_vel - V_REFF)**2

        W4 = np.sign(-v_rel) * MS * 9.81 * ALPHA

        W = W1 + W2 + W3 + W4

        thrust_Factor_a = 16.26;
        thrust_Factor_b = 6.25;
        thrust_Factor_c = 0.8;

        # Compute thrust of vessel
        thrust = thrust_Factor_c * (thrust_Factor_a * self.ships.power[self.id]) / (self.ships.power[self.id]**(1/3) + thrust_Factor_b * v_rel)

        a_next = 1000
        a_prev = 1000

        # Compute acceleration
        dv_str_dt = (self.obs_stream_vel - prev_stream_vel) / TIMESTEP

        mh_prev = self.hydrodyn_mass(a_prev, AS, MS)
        mh_next= self.hydrodyn_mass(a_next, AS, MS)
        dmh_dt = (mh_next - mh_prev) / TIMESTEP

        # Bewegungsgleichung (Motion equation)
        new_acc = 1.0 / ((MS + MH) * (thrust - W - v_rel * dmh_dt + MH * dv_str_dt))
        self.wid = new_acc

        # Calculate CF Parameter
        CF_GAMMA = 0.82
        CF_GAMMA_L = 0.8
        C_AST = 3.2
        C_I = 0.5
        F_V2 = 0.3
        C_KORR = 1.4
        BETA = 2.0
        BETA_C = 0.5
        DELTA_C_FS = 0
        DELA_C_FS_KORR_BS = 0

        cd_0 = 0.22 * (LS / (CF_GAMMA * BS))**BETA_C
        cd_max = cd_0 + (C_AST -cd_0) * (TSAVG / h)**BETA
        vStr_vSdW = self.obs_stream_vel / v_rel
        F1 = (1.0 + 2.0 * vStr_vSdW) + F_V2 * (vStr_vSdW)**2

        F11 = (TSAVG / h) / (1.0 - (TSAVG / h))
        F2 = 1.0 + 0.35 * np.log(BS / TSAVG + 1.0)
        F3 = 0.5 * np.pi * TSAVG**2 / (TSAVG * BS + 0.2 * TSAVG * TSAVG)
        F4 = (TSAVG * BS + 0.2 * TSAVG**2) / (TSAVG * BS)
        CMHY = F4 * (F3 * F2 + F11)

        k_prime = 1.0 - (3.4^(-0.22 * (LS / BS - 1.0)))
        CF = 0.5 + C_I * (BS / LS) + (2.0 / cd_max * BS / LS * CF_GAMMA / CF_GAMMA_L * C_KORR * F1 + 1.0 / 6.0) /\
             (CF_GAMMA / CF_GAMMA_L * BS / LS * CMHY / cd_max * 4 * k_prime + 1.0)

        CF += DELTA_C_FS + DELA_C_FS_KORR_BS

        # This returns an immutable tuple as it is faster
        return (new_acc, new_squat, CF)

    @staticmethod
    def hydrodyn_mass(A,AS,MS):
        
        N = A/AS
        MASS = MS * (1/(N - 1) + 0.1)

        return MASS

river = River()
ships = Ships(river=river, num_ships=5, ship_lengths= [100]*5, ship_widths= [10]*5,ship_mass=[1e5]*5,y_location=[150]*5)
lat = LatConPol(0)
f = lat.compute_obs(ships, river,2)