import numpy as np
import os
from IPython.display import HTML
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from constants import *
from states_functions import *

def load_points(filename):
    """Loads the 2 columns from the path \Bachelor\Steady States\\filename and returns them as a 2D numpy array."""
    current_folder = os.getcwd()
    parent_folder = os.path.dirname(current_folder)
    sub_folder = "\\Steady States\\Cell States\\"
    file_path = parent_folder + sub_folder + filename
    
    # Load the file using numpy
    data = np.loadtxt(file_path, usecols=(0, 1))

    # Access the columns
    column1 = data[:, 0]
    column2 = data[:, 1]
    points = np.column_stack((column1, column2))
    return points


def save_protein_steady_state(proteins, protein_names, filename):
    """Saves the steady state of the proteins in a text file"""
    with open(filename, "w") as file:
        
        # Write the protein names
        for i, name in enumerate(protein_names):
            file.write(name + " ; ")
            # Writes the last frame of each proteins value as a single long list:
            for j, level in enumerate(proteins[i][-1]):
                if j == len(proteins[i][-1])-1:
                    file.write(str(level))
                else:
                    file.write(str(level) + ",")
                
            file.write("\n")
    file.close()
    return None

def load_protein_steady_state(filename):
    """Loads the steady state of the proteins from a text file"""
    current_folder = os.getcwd()
    parent_folder = os.path.dirname(current_folder)
    sub_folder = "\\Steady States\\Protein Levels\\"
    file_path = parent_folder + sub_folder + filename
    
    with open(file_path, "r") as file:
        lines = file.readlines()
        proteins = []
        for line in lines:
            protein = line.split(" ; ")[1]
            protein = protein.split(",")
            protein_list = []
            for level in protein:
                protein_list.append(float(level))
            proteins.append(protein_list)
    return np.array(proteins[0]), np.array(proteins[1])


def get_steady_state(points_from_euler, filename = "SteadyState.txt"):
    """After en Euler simulation we should save the locations for
    the cells so we don't need to keep simulating their movement but can load a file with points"""
    points = points_from_euler[-1] 
    file = open(filename, "w") 
    file.write("Now the file has more content!") 
    file.close()
    np.savetxt(filename, points) # Save the points to a file


def gen_points(N_points, R = 1):
    """Generates N randomly placed "Cells" in a numpy array"""
    points = []
    # Generates random points
    for _ in range(N_points):
        r = R * np.sqrt(np.random.rand())
        theta = np.random.rand() * 2 * np.pi
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        location = np.array([x,y])
        points.append(location)
    return np.array(points)


def generate_spiral_points(N, rad_increment=0.1, angle_increment=0.1):
    """Generates N points in a spiral pattern and returns them as a 2D numpy array."""
    points = np.zeros((N, 2))
    angle = 0
    radius = 0.1
    
    for i in range(N):
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        points[i] = [x, y]
        
        angle += angle_increment
        radius += rad_increment
    
    return points


### Animation functions ###


def update(i, list_of_SMADS, list_of_NOGGINS, t_list, points, fig, axis, vor):
    """
    Update function is what is happening every frame. Imagine that you have to make a new plot every frame,
    but instead you clear the background and overlay something else.
    """
    # Step 1 # Generate data
    SMAD = list_of_SMADS[i]
    NOGGIN = list_of_NOGGINS[i]
    t = t_list[i]
    
    
    # Step 2 # Clear previous plot and plot the updated plot
    axis.clear()          # Alternatively ax.clear() Clears the entire thing, but collections offers selection

    for i, point in enumerate(points):
        circle_color = color_check(i, colors, SMAD, NOGGIN)
        
        # Plots the circle
        circle = plt.Circle((point[0], point[1]), r_cell, facecolor = circle_color, edgecolor = "black", alpha = 0.8)
        axis.add_patch(circle)

    # Step 3 # Modification of axis layout should occur in here
    axis.axis("equal")
    axis.set_xlim(screen_size)
    axis.set_ylim(screen_size)
    return None



def color_update(i, list_of_SMAD_levels, list_of_NOGGIN_levels,
                 points, axis, vor, colormap, norm,
                 colormap_mode = "value", protein_to_plot = "SMAD"):
    """
    Update function is what is happening every frame.
    Imagine that you have to make a new plot every frame,
    but instead you clear the background and overlay something else.
    
    Arguments:
    axis: The axis to plot on
    vor: The voronoi class
    colormap: The colormap to use (Initialized beforehand)
    norm: The normalization of the colormap (Initialized beforehand)
    colormap_mode: The mode of the colormap, either "value", "mean"
    protein_to_plot: The protein to plot, either "SMAD" or "NOGGIN" or later "WNT" and "NODAL"
    """
    # Step 1: Generate data
    SMAD = list_of_SMAD_levels[i]
    NOGGIN = list_of_NOGGIN_levels[i]
    
    # Step 2: Clear previous plot and plot the updated plot
    axis.clear()  # Alternatively ax.clear() Clears the entire thing, but collections offers selection


    # Protein selection
    if protein_to_plot == "SMAD":
        protein = SMAD
    elif protein_to_plot == "NOGGIN":
        protein = NOGGIN
    else:
        raise ValueError("Protein not recognized")

    # Mode selection
    if colormap_mode == "value": 
        for i, point in enumerate(points):
            # Now using the SMAD levels we need to generate a color intensity for the circle
            circle_color = colormap(norm(protein[i]))  # Map SMAD levels to colors using the colormap
            
            # Plot the circle
            circle = plt.Circle((point[0], point[1]), r_cell, 
                                facecolor=circle_color, edgecolor="black", alpha=0.8)
            axis.add_patch(circle)
    elif colormap_mode == "mean":
        for i, point in enumerate(points):
            # Generate the mean levels
            SMAD_level, SMAD_mean = get_level_of_protein(protein, i, vor)
            
            # Now using the SMAD levels we need to generate a color intensity for the circle
            circle_color = colormap(norm(protein[i]))  # Map SMAD levels to colors using the colormap
            
            # Plot the circle
            circle = plt.Circle((point[0], point[1]), r_cell,
                                facecolor=circle_color, edgecolor="black", alpha=0.8)
            axis.add_patch(circle)

    # Step 5: Modification of axis layout should occur here
    axis.axis("equal")
    axis.set_xlim(screen_size)
    axis.set_ylim(screen_size)
    return None





def animate(fig, update_func, func_arguments, n_frames, interval = 300,
            fps = 4, save_anim = False, gif_name = "Test.gif"):
    """Creates a basis for animating so it doesn't fill as much.
    We can add different update functions and also initialize
    colormaps prior to initializing this function
    
    The figure and other stuff that isn't being updated needs to be initialized prior and passed as arguments
    """
    
    # Running the animation
    anim = FuncAnimation(fig, update_func, fargs=func_arguments,
                         frames=n_frames, interval=interval) # Controls the plotted animation
    plt.rcParams['animation.embed_limit'] = 2**128          # Limits the Bytes allocated to the animation
    if save_anim:
        anim.save(gif_name, writer='Pillow', fps=fps)      # Controls the animation for the downloadable gif and saves it
    plt.close(fig)
    return HTML(anim.to_jshtml())      # Converts anim to an actual gif

def color_check(cell_index, colors, SMAD, NOGGIN):
    """ A very bad function for picking out colors dependent on states
    """
    is_SMAD = SMAD[cell_index] == 1
    is_NOGGIN = NOGGIN[cell_index] == 1
    if (is_SMAD and is_NOGGIN):
        color = colors[2]

    elif (is_SMAD  and not is_NOGGIN):
        color = colors[1]

    elif (not is_SMAD and is_NOGGIN):
        color = colors[3]

    else:
        color = colors[0]
        
    return color