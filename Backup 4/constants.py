import numpy as np
from my_functions import *

# Initial conditions
S = 1
beta = 5
dt = 0.5
t_start = 0
t_end = 300
n_points = 625
screen_size = (-30, 30)


# Threshhold values
theta_SMAD = 0
theta_NOGGIN = 7
theta_WNT = 0
theta_NODAL = 0

threshholds = [theta_SMAD, theta_NOGGIN, theta_WNT, theta_NODAL]


# Cell states:
BMP = 1     # BMP is everywhere in the cell
# Creating cell lists
SMAD = np.zeros(n_points)
NOGGIN = np.zeros(n_points)


# Changing the noise
mu = 0
sigma = 0.025

# Getting the size of the cell (Plotting)
r_cell = get_cell_radius(beta, S)

colors = ["blue", "red", "yellow", "green"]     # We try and differentiate the different cells 
                                                # blue = NONE, red = SMAD+, yellow = SMAD,NOGGIN , green = NOGGIN
