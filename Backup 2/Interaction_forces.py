# A file that contains the functions for the interaction forces between the cells
import numpy as np
from my_functions import *
import math as m
from constants import *
from Voronoi_functions import *
from scipy.spatial import Voronoi, voronoi_plot_2d

def gaussion_noise(N, mean=mu, std=sigma, dimensions = 2):
    """Generates random noise in the size of N
    
    Args:  
        N (int): the size of the noise
        mean (float): the mean of the noise
        std (float): the standard deviation of the noise ; a number for how far the cells can move
        dimensions (int): (x,y) -> 2, (x,y,z) -> 3 (dependent on the number of dimensions in the problem)
    
    """
    return np.random.normal(mean, std, size = (N,dimensions))


def interaction_force(point1, point2, beta = 5, S = 1):
    """Generates the interaction force between two cells and splits up into carteisan coordinates

    Returns:
        list: the components of the force in (F_x, F_y)
    """
    dist = m.dist(point1, point2)
    r_vec = (point1 - point2)/dist
    force = np.exp(-dist) - (S/beta)*np.exp(-dist/beta)
    
    force_x, force_y = force * r_vec
    return np.array([force_x, force_y])


def get_cell_forces(points, beta=5, S=1):
    """
    Returns:
        list of list: Returns a combined nested list where first element is all the forces components on the first cell.
        Each sublist is composed of (F_x, F_y) 
    """
    
    # The idea is that we should generate the neighbours before we calculate the forces
    # This means that we only calculate the forces between the cells that are neighbours
    # To do this we shuld generate the Voronoi diagram and obtain the neighbours
    # Then find a nice way to loop through the interactions
    
    # Voronoi
    vor = Voronoi(points)
    
    # We start by only checking the neighbours and disregard the nearest neighbours for now
    neighbours = get_nearest_neighbours(vor)
    
    forces = []                         # All forces for all points
    
    # For every point we calculate the forces from the neighbours on that point
    # and returns it as a single value, which is done for all points
    for i, point1 in enumerate(points):
        cell_forces = []                       # The forces for each cell
        for index in neighbours[i]:            # Gets the neighbours index for the specific cell                  # For every neighbour belonging to that cell                                     # Not counting the interaction of the cell on itself (even though its zero)
            force = interaction_force(point1, points[index], beta, S)
            cell_forces.append(force)
        
        # The idea is that we sum all the x_components together to form one interaction force
        force_on_cell = np.sum(cell_forces, axis=0)         # [sum(x_1, x_2, ..., x_n) and sum(y_1, y_2, ..., y_n)]
        forces.append(force_on_cell)
    return np.array(forces)