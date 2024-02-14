import numpy as np
import os

def load_points(filename):
    """Loads the 2 columns from a file and returns them as a 2D numpy array."""
    current_folder = os.getcwd() + "/Steady States/"
    file_path = current_folder + filename
    
    # Load the file using numpy
    data = np.loadtxt(file_path, usecols=(0, 1))

    # Access the columns
    column1 = data[:, 0]
    column2 = data[:, 1]
    points = np.column_stack((column1, column2))
    return points


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


# With this we should for each point generate a radius
def get_cell_radius(beta = 5, S = 1):
    """At equillibrium the distance between two cells are equal to the following 
    which then needs to be divided by 2 to get one cell"""
    return (beta*np.log(beta/S))/(beta-1) * (1/2)


def color_check(cell_index, colors, SMAD, NOGGIN):
    """ A very bad function for picking out colors dependent on states
    """
    is_SMAD = SMAD[cell_index]
    is_NOGGIN = NOGGIN[cell_index]
    if (is_SMAD and is_NOGGIN):
        color = colors[2]

    elif (is_SMAD and not is_NOGGIN):
        color = colors[1]

    elif (not is_SMAD and is_NOGGIN):
        color = colors[3]

    else:
        color = colors[0]
        
    return color