# Molecular dynamics (MD) and Langevin integrators for MolSSI summer school.

# 2017-7
# Basic setting

import numpy as np
import scipy.integrate as sc
from scipy.stats import norm
import copy
import math

# Parameters
r_sw = 9.0
r_cut = 10.0
timestep = 0.0001
gammma = 1/(500*timestep)   #gamma is the collision rate
kbT = 1
mass = 1.0
particle_num = 10
box_length = 10.0

# import coodinates
coordinate = np.loadtxt("lj_sample_config_periodic1.txt", skiprows=2,usecols=(1,2,3))
velocity = np.zeros(particle_num*3).reshape(particle_num,3)

# MD simulation

total_step = 10
total_energy = 0.0
mass = 1.0

frc = np.arange(3*particle_num).reshape(particle_num,3)

file = open( "output.x", "w" )
file.close()
file_2 = open ("output.v","w")
file.close()

file = open( "output.x", "ab" )
file_2 = open ("output.v","ab")

initial_velocities = []
#initialize velocity
for atom in range(0, particle_num):
    v_initial = initialize_velocity(mass)
    initial_velocities.append(v_initial) 
    
for step in range(0, total_step):
    
    np.savetxt(file, coordinate, fmt='%10.6f')
    np.savetxt(file_2, velocity, fmt='%10.6f')
    
    coordinate_new = np.zeros(particle_num*3).reshape(particle_num,3)
    velocity_new = np.zeros(particle_num*3).reshape(particle_num,3)
    
    for  i in range(0, particle_num):
        
        #OVRVO algorithom
        x = coordinate[i].copy()
        v = initial_velocities[i]
        v_1 = np.exp(- 0.5*gammma*timestep)*v + np.sqrt(1-np.exp(-gammma*timestep))*np.sqrt(kbT/mass)*np.random.rand()
        v_2 = v_1 + 0.5*timestep/mass*force(i)
        x_new = x + timestep*v_2
        coordinate_new[i] =  x_new
        v_3 = v_2 + 0.5*timestep/mass*force(i)
        v_new = np.exp(- 0.5*gammma*timestep)*v_3 + np.sqrt(1-np.exp(-gammma*timestep))*np.sqrt(kbT/mass)*np.random.rand()
        velocity_new[i] =  v_new
    
    coordinate = coordinate_new.copy()
    velocity =  velocity_new.copy()
    
    print(str(total_step)+ ' steps in total, finished ' + str(step+1) + ' steps')
        
file.close()

file_2.close()
        
