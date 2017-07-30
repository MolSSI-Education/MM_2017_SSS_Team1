# Molecular dynamics (MD) and Langevin integrators for MolSSI summer school.

# 2017-7
# Basic setting

import numpy as np
import scipy.integrate as sc
from scipy.stats import norm
import copy
import math

#FUNCTIONS
# L-J Potential & the scaler cooeficient of its gradient

def lennard_jones_potential(r):
    sig_by_r6 = np.power(1/r,6)
    sig_by_r12 = np.power(sig_by_r6,2)
    return 4.0*(sig_by_r12-sig_by_r6)

# it needs to be multipled by the vector r_ki
def lj_gradient_co(r):
    sig_by_r8 = np.power(1/r,8)
    sig_by_r14 = np.power(1/r,14)
    return -48.0*(sig_by_r14 - 0.5*sig_by_r8)

# Switching function & the scaler cooeficient of its gradient
def switching_function(r):
    if r <= r_sw:
        return 1.0
    
    elif r < r_cut:
        r_2 = np.square(r)
        r_cut_2 = np.square(r_cut)
        r_sw_2 = np.square(r)
        s = np.square((r_cut_2-r_2))*(r_cut_2+2*r_2-3*r_sw_2)/np.power((r_cut_2-r_sw_2),3)
        return s
    
    else: return 0.0

# it needs to be multipled by the vector r_ki
def s_gradient_co(r):
    
    if r <= r_sw:
        return 0.0    
    elif r < r_cut:
        r_2 = np.square(r)
        r_cut_2 = np.square(r_cut)
        r_sw_2 = np.square(r)
        co = 12.0*(r_2-r_cut_2)*(r_2-r_sw_2)/np.power((r_cut_2-r_sw_2),3)
        return co
    else: return 0.0
def single_force(r_vec, r):
    return - (switching_function(r)*lj_gradient_co(r) 
              + lennard_jones_potential(r)*s_gradient_co(r))*r_vec + Long_range_dispersion(r)
    
    
    
def force(i_particle):
    frc = np.zeros(3)   
    for j in range(0, particle_num):
        if j != i_particle:
            vec_rij = coordinate[j] - coordinate[i_particle]
            vec_rij = vec_rij - box_length*np.round(vec_rij/box_length)
            #periodic boundary
            scal_rij = np.linalg.norm(vec_rij)
            frc += single_force(vec_rij, scal_rij) 
    
    return frc

def integrand(r):
     return  r ** 2 

def initialize_velocity(m):
    variance = kbT / m
    normal1 = norm(0, variance)
    normal2 = math.sqrt(1 / 2 * 3.1415 * variance)
    mu = 0
    sampled_velocity = np.random.normal(mu, variance)
    return sampled_velocity

def Long_range_dispersion(rij):
    rij2 = rij * rij
    U_lj = lennard_jones_potential(rij)
    U_tail = tail_correction(box_length)
    S_ = switching_function(rij)
    integral1, integral2 = sc.quad(integrand, r_sw, r_cut)
    integral = integral2 - integral1
    lrj = U_tail + integral * U_lj * (1-S_)
    return lrj

def tail_correction(box_length):
    volume = np.power(box_length,3)
    sig_by_cutoff3 = np.power(1.0/r_cut, 3)
    sig_by_cutoff9 = np.power(sig_by_cutoff3, 3)
    e_correction = sig_by_cutoff9 - 3.0 * sig_by_cutoff3
    e_correction *= 8.0 / 9.0 * np.pi * particle_num / volume * particle_num
    return e_correction



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
        
