from constants import *
from my_functions import *
from neighbour_functions import *
from states_functions import *

def diffusion_step(SMAD, Voronoi, n_steps = 10, select_noSMAD_cells = False):
    """A function that should be used in collaboration with an euler loop. 
    Currently the diffusion is done by choosing an index of a cell that
    is giving its protein to its neighbours. This means that the proteins
    maybe are diffusing rather fast in the beginning which might not be accurate
    
    Returns:
        SMAD (list): An updated list of protein values
    """
    
    # We repeat the step n_steps times
    for _ in range(n_steps):
        # We choose a random cell that has a non-zero protein value
        possible_cells = np.where(SMAD > 0)[0]         # np.where picks out the indices that satisfies the condition given
                                                       # Don't really need to use np.where, but it is a good practice # Keep this in mind
                                                       # np.where returns a tuple so what we need is in the first index
        
        if select_noSMAD_cells:
            ### I Feel like Code is horribly slow ###        
            # A for loop that only picks out cells with more than 3 neighbours so we don't
            # select diffusion that occurs Where diffusion most likely wouldn't change anything
            neighbours_of_possible_cells = [get_nearest_neighbours(Voronoi.ridge_points, cell, Voronoi.points) for cell in possible_cells]
            
            new_possible_cells = []
            for i, neighbours in enumerate(neighbours_of_possible_cells):   # We check all cells and see if their neighbours are only SMAD
                onlySMAD = True
                for neighbour in neighbours:                                # If they are, we don't want to diffuse from that cell
                    if SMAD[neighbour] < 1:
                        onlySMAD = False
                # If the cell has neighbours that are not only SMAD, we add it to the new list        
                if not onlySMAD:
                    new_possible_cells.append(possible_cells[i])
                    
            possible_cells = new_possible_cells
    
    
    
        if len(possible_cells) == 0:                # We check that there are cells that can diffuse
            raise ValueError("No cells with protein to diffuse, List empty!")
    
        # Then pick one of these cells to diffuse
        cell_index = np.random.choice(possible_cells)

        # We get the cells neighbours
        neighbours = get_nearest_neighbours(Voronoi.ridge_points, cell_index, Voronoi.points)
        
        # We can only diffuse to cells with lower values than ourself
        diffusive_neighbours = np.where(SMAD[neighbours] < SMAD[cell_index])[0] # This picks out the index in the neighbours list that fulfills the condition
        
        # What should diffuse to what neighbours
        protein_to_remove = 0
        for index in diffusive_neighbours:
            neighbour = neighbours[index]
            
            # We calculate the fraction of the protein that should be given
            protein_to_give = ((SMAD[cell_index] - SMAD[neighbour])/SMAD[cell_index]) * SMAD[cell_index]/(len(diffusive_neighbours)+1)
            
            # The reasoning for above is that we want a diffusion whose strength is dependent
            # on how big the concentration of the neighbours are. Therefore we take a difference to get a strength factor
            # and then divide by the number of neighbours we can diffuse to so we split up the protein and don't kill the cell
            
            # We add the protein to neighbours and remove it from the cell
            SMAD[neighbour] += protein_to_give
            protein_to_remove += protein_to_give
            
        SMAD[cell_index] -= protein_to_remove
    return SMAD
