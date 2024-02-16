import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial import Voronoi
from neighbour_functions import *
from constants import *

# This function below should be used in collaboration with an euler loop to get the changing cell states
# each timestep. This should then be saved and be made into an animation.
def state_of_cell(cell_index, Voronoi, SMAD, NOGGIN, threshholds):
    """A function that changes the state of the cell when certain conditions are met.
    Some of these can be changing from SMAD- to SMAD+ which is saved as a change in the
    SMAD list: 
    SMAD[1] = 0 ; The second cell has SMAD-
    SMAD[1] = 1 ; The second cell has SMAD+

    Args:
        cell_index (int): The index of the cell we are looking at
        Voronoi (Voronoi): The created Voronoi class that holds points, neighbours and so on
        SMAD (list): List of cells that are SMAD
        NOGGIN (list): List of cells that are NOGGIN
        threshholds (list): list of the different threshhold values (SMAD, NOGGIN, WNT, NODAL)
    """
    threshhold_SMAD = threshholds[0]
    threshhold_NOGGIN = threshholds[1]
    
    # We first need to get the indicies of the neighbours of the cell we are looking at to determine the state
    nearest_neighbours = get_nearest_neighbours(Voronoi.ridge_points, cell_index, Voronoi.points)
    
    # Calculate the levels of the different proteins in the cell and its neighbours
    NOGGIN_levels = np.sum(NOGGIN[nearest_neighbours]) + NOGGIN[cell_index]     # This indexing works because
    SMAD_levels = np.sum(SMAD[nearest_neighbours]) + SMAD[cell_index]           # they are numpy arrays
    
    # Convert to SMAD
    if BMP * (len(nearest_neighbours)+1 - NOGGIN_levels) > threshhold_SMAD:
        SMAD[cell_index] = 1        # Converts to SMAD+
    else:
        SMAD[cell_index] = 0        # Converts to SMAD-
    
    
    # Convert to NOGGIN
    if (NOGGIN_levels + SMAD_levels) > threshhold_NOGGIN:
        NOGGIN[cell_index] = 1      # Converts to NOGGIN+
    
    else:
        # Generate a random number and convert NOGGIN to NOGGIN-
        prob = np.random.uniform(0,1) 
        if prob < 0.25:
            NOGGIN[cell_index] = 0
        else:
            NOGGIN[cell_index] = NOGGIN[cell_index]
    
    return SMAD, NOGGIN

def get_levels_of_protein(protein, cell_index, voronoi_class):
    neighbours = get_nearest_neighbours(voronoi_class.ridge_points, cell_index, voronoi_class.points)
    levels = np.sum(protein[neighbours]) + protein[cell_index]
    return levels