import numpy as np
from scipy.spatial import Voronoi

### Voronoi Functions ###

### VORONOI CLASS ###
# this class contains the following attributes:
# - points: the input points (de punkter som jeg selv har givet)
# - vertices: Krydsningerne mellem linjerne i Voronoi plottet
# - ridge_points: en liste af punkter over de celler der deler en linje (alle interaktioner)

def check_line_of_sight(cell_index, neighbour_index, local_neighbours_index, points):
    """Finds the center between two points and checks in any of the neighbours in the local neighbours
       are closer to the center. 
       If Yes then the interaction is not a nearest neighbour.

    Returns:
        bool: True if the interaction is a nearest neighbour, False if not
    """
    # Calculate center point
    center = (points[cell_index] + points[neighbour_index])/2
    distance_from_cell = np.linalg.norm(center - points[cell_index])
    
    # Get distance from center to some of the other points
    # Vi looper igennem alle interaktioner med det indeks vi har fået og tjekker om nogle af de punkter
    # ligger tættere og dermed forstyrrer eye of sight af de andre punkter
    check = True
    for neighbour in local_neighbours_index:
        if neighbour == neighbour_index:
            break
        # Calcs distances between the two points
        distance_from_neighbour = np.linalg.norm(center - points[neighbour])

        # if a < b then the interaction is a nearest neighbour
        if distance_from_cell > distance_from_neighbour:
            check = False
    
    return check


# Github kode for at udvælge naboer (FAST)
def get_interacting_neighbors(ridge_points, interaction_index):   
    """Checks all the ridge points in the Voronoi class and sorts the interacting neighbours for a given cell cell_index
       Interacting just means that the two cells share a line

    Args:
        ridge_points (list): list of ridgepoints (should be in the form of Voronoi.ridge_points)
        interaction_index (int): The cell_index of the cell we want to find the interacting neighbours for

    Returns:
        list: A list of all the interacting cell indices for a given cell cell_index
    """
    # Convert ridge_points to NumPy array for faster computations
    ridge_points = np.array(ridge_points)
    
    # Find interactions involving the given point cell_index
    # Looper igennem alle punkter og tjekker om indexet vi fokuserer på er i listen som det første element eller det andet
    interactions = ridge_points[np.logical_or(ridge_points[:, 0] == interaction_index, ridge_points[:, 1] == interaction_index)]
    
    # Fjerner indexet fra interaktionerne så vi kun har naboerne i en liste
    # [[1,2], [1,3], [1,4]] -> [2,3,4]
    interactions = interactions[interactions != interaction_index]
    
    return np.array(interactions)


def get_nearest_neighbours(voronoi_interactions, cell_index, points):
    """Gets the neighbours of all points and checks if the neighbours are classified as "nearest neighbours"
    according with Maja's definition. (b > a).
    By Nearest we mean that the two cells are closer to each other than to any other cell."""
    # Loops through all the points, gets the corresponding neighbours and checks if those neighbours are considered nearest
    nearest_neighbours = []
    local_neighbours_index = get_interacting_neighbors(voronoi_interactions, cell_index)

# Checks if the interaction is a nearest neighbour
    for neighbour in local_neighbours_index:
        is_nearest_neighbour = check_line_of_sight(cell_index, neighbour, local_neighbours_index, points)
    
        # Adds the nearest neighbours to the list
        if is_nearest_neighbour:
            nearest_neighbours.append(neighbour)

    # Eventually we need to return more than just the neighbours (maybe points)
    return nearest_neighbours

