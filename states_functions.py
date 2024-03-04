from scipy.spatial import Voronoi
from neighbour_functions import *
from constants import *
from my_functions import *
from diffusion import diffusion_step
import copy 

def reduce_state(protein_states, protein_levels, neighbours, cell, threshhold = 1):
    """This function produces two things. Firstly it takes the sum of the protein levels and sees if it is able
    to convert the cell to a different state. It also does this for the cell index that we are given"""
    total_neighbours = protein_states[neighbours] + protein_levels[neighbours]
    indices = np.where(total_neighbours > threshhold)
    total_neighbours[indices] = threshhold
    
    total_cell = protein_states[cell] + protein_levels[cell]
    if total_cell > threshhold:
        total_cell = threshhold
    
    return total_neighbours, total_cell


# This function below should be used in collaboration with an euler loop to get the changing cell states
# each timestep. This should then be saved and be made into an animation.
def state_of_cell(cell_index, Voronoi, SMAD, NOGGIN, SMAD_levels, NOGGIN_levels, threshholds):
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
    
    # We make an intermediate step to see if concentration is high enough for convertion:
    SMAD_neighbours, SMAD_cell = reduce_state(SMAD, SMAD_levels, nearest_neighbours, cell_index)
    
    # Calculate the levels of the different proteins in the cell and its neighbours
    NOGGIN_amount = np.sum(NOGGIN[nearest_neighbours]) + NOGGIN[cell_index]     # This indexing works because
    SMAD_amount = np.sum(SMAD_neighbours) + SMAD_cell           # they are numpy arrays
    SMAD_test = np.sum(SMAD[nearest_neighbours]) + SMAD[cell_index]
    
    # Convert to SMAD
    if BMP * (len(nearest_neighbours)+1 - NOGGIN_amount) > threshhold_SMAD:
        SMAD[cell_index] = 1        # Converts to SMAD+
    else:
        SMAD[cell_index] = 0        # Converts to SMAD-
    
    
    # Convert to NOGGIN
    if (NOGGIN_amount + SMAD_amount) > threshhold_NOGGIN:
        NOGGIN[cell_index] = 1      # Converts to NOGGIN+
    
    else:
        # Generate a random number and convert NOGGIN to NOGGIN-
        prob = np.random.uniform(0,1) 
        if prob < noggin_conversion:
            NOGGIN[cell_index] = 0
        else:
            NOGGIN[cell_index] = NOGGIN[cell_index]
    
    return SMAD, NOGGIN

def test_update_protein_levels(protein, protein_levels, cell):
    """A test function for checking if Montecarlo simulation should be used for the protein levels"""
    protein_levels[cell] = protein_levels[cell] + protein[cell]
    if protein_levels[cell] > 1:
        protein_levels[cell] = 1
    return protein_levels

def update_protein_levels(protein, protein_levels):
    """We should make a function that updates the protein levels of the cells. 
    For now this function updates all cells"""
    # We imagine that the cell is generating a set amount of protein each timestep
    protein_levels =+ protein
    # if it is above one, we set it to one
    indices = np.where(protein_levels > 1)
    protein_levels[indices] = 1
    
    return protein_levels
    
def diffusion_boundary_check(protein_levels):
    """Checks that the protein levels aren't becoming negative"""
    indices = np.where(protein_levels < 0)
    protein_levels[indices] = 0
    return protein_levels



def get_level_of_protein(protein, cell_index, voronoi_class):
    """Gets the levels of a specific protein in a cell and returns it as a sum

    Args:
        protein (list): A list of the protein levels in the cells, should be the level and not the lists holding the states
        cell_index (int): Trivial
        voronoi_class (Voronoi): The class that holds the voronoi information

    Returns:
        level, mean: The sum of the protein in the cell and the average level of the protein in the cell
    """
    # We generate the list of nearest neighbours
    neighbours = get_nearest_neighbours(voronoi_class.ridge_points, cell_index, voronoi_class.points)
    # Obtain the amount of the given protein in the cell and its neighbours
    level = np.sum(protein[neighbours]) + protein[cell_index]
    
    # Calculate the mean level and remembering that the cell itself is also a neighbour
    mean = level/(len(neighbours) + 1)
    return level, mean



def simulate_BMP(t_start, t_end, dt, points, SMAD, NOGGIN, SMAD_levels, NOGGIN_levels, self_degredation = 0.01, diffusion = True, diffusive_steps = 1000, diffusion_mode = "simple"):
    """Simulates the BMP interactions between cells.
    It is assumed that BMP is supplied to the cells during the whole simulation. 
    The simulation should only look at the interactions but not move the cells

    Arguments:
    points (list): List of Steady state points
    diffusion (bool): If the simulation should include diffusion
    diffusive_steps (int): How many times the diffusion should occur each timestep
    diffusion_mode (str): The mode of diffusion. "simple" or "noSMAD"


    Returns:
        lists: Yes
    """
    # Obtaining initial vals
    vor = Voronoi(points)
    SMAD_lists = [copy.deepcopy(SMAD)]
    NOGGIN_lists = [copy.deepcopy(NOGGIN)]
    SMAD_levels_lists = [copy.deepcopy(SMAD_levels)]
    t_list = [t_start]
    
    # Running the euler simulation
    while t_start < t_end:
        
        # Calculate changes in states and values for each timestep
        for _ in range(len(points)):    # Or For each point we check a random cells state and update it
            cell = np.random.randint(0, len(points)) # Generates index from 0 to N-1
            
            # Here we update the state of the index before letting the 
            # other cells update according to the previous timestep? Is this a Good IDEA?!
            SMAD, NOGGIN = state_of_cell(cell, vor, SMAD, NOGGIN, SMAD_levels, NOGGIN_levels, threshholds)
            
        
        # Update the protein levels
        SMAD_levels = update_protein_levels(SMAD, SMAD_levels)
        
        # We try and do another monte carlo simulation but off the protein levels
        # for _ in range(len(points)):
        #     cell = np.random.randint(0, len(points))
        #     SMAD_levels = test_update_protein_levels(SMAD, SMAD_levels, cell)
        
        
        # We let the system diffuse to equillibrium
        if diffusion:
            SMAD_levels = diffusion_step(SMAD_levels, vor, n_steps = diffusive_steps, select_noSMAD_cells = False)
                
        # Lastly we remove a set amount from every cell of each protein because of self degradation
        SMAD_levels -= self_degredation
        SMAD_levels = diffusion_boundary_check(SMAD_levels)
        
        # Update time step
        t_start += dt
        
        # Update lists
        SMAD_lists.append(copy.deepcopy(SMAD))
        NOGGIN_lists.append(copy.deepcopy(NOGGIN)) 
        SMAD_levels_lists.append(copy.deepcopy(SMAD_levels))
        t_list.append(t_start)  
    
    return t_list, SMAD_lists, NOGGIN_lists, SMAD_levels_lists