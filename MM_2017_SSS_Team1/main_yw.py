import numpy as np
import scipy as sc

from mpl_toolkits.mplot3d import Axes3D
%matplotlib notebook
coordinates = np.loadtxt("lj_sample_config_periodic1.txt", skiprows=2, usecols=(1,2,3))
num_particles = len(coordinates)
initial_velocities = {}
# #initialize velocities
for atom in range(0, num_particles):
    v_initial = initial_velocity() 
    initial_velocities[coordinates[atom]] = v_initial
    print(initial_velocities)
#MD Algorithm 

coordinates_of_simulation = []
    for atom in range (0, num_particles):
        print("Simulation for particle", atom + 1)
        integrator_algorithm(num_steps,v_initial, coordinates[atom]) #Grey. Return X_t_plus_1 and current velocity to be used for calculation of next step. and put coordinates in array, should return coordinates. of simulationn. 
 
    
