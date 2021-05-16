"""
A ForceField for only Ar molecules. Due no chemical bonding, only non-bonding terms are working out. Inspired in the AMBER ForceField, it will use the Lennard-Jones 12-6 for Van der Waals interactions, and columbic interactions.
Ar-Ar parameters -> https://courses.physics.illinois.edu/phys466/sp2013/projects/2001/team5/node13.html
"""
import numpy as np

epsilon = 119.8*1.380649*(10**-23) # Ar-Ar value
sigma = 3.41 #Ar-Ar value
charge = 0.5 #Atom charge. Must be changed because it is an aproximation
varepsilon = 1 * 8.8542*10**-12 #Permeability of medium. If first multiplication term is 1, it means we're in vacuum

def PotentialEnergy(positions):
    global epsilon, sigma, charge, varepsilon
    _Nparticles = len(positions)
    V = 0
    for i in range(0, _Nparticles-1):
        for j in range(i+1, _Nparticles):
            r = np.sqrt((positions[str(i)][0]-positions[str(j)][0])**2 + (positions[str(i)][1]-positions[str(j)][1])**2 + (positions[str(i)][2]-positions[str(j)][2])**2)
            V += (4*epsilon*( (sigma/r)**12 - (sigma/r)**6)) + ((charge**2) / (varepsilon*r))       
    return V

def Force(positions, candidate):
    #Force that actues over a candidate particule
    global epsilon, sigma, charge, varepsilon
    _Nparticles = len(positions)
    Fx, Fy, Fz = 0, 0, 0
    for i in range(0, _Nparticles):
        if i == candidate:
            continue
        else:
            rx = np.sqrt((positions[str(i)][0]-positions[str(candidate)][0])**2)
            ry = np.sqrt((positions[str(i)][1]-positions[str(candidate)][1])**2)
            rz = np.sqrt((positions[str(i)][2]-positions[str(candidate)][2])**2)
            Fx += (4*epsilon* ((12*((sigma**12)/(rx**13))) - (6*( (sigma**6) / (rx**7))) )) + ((charge**2) / (varepsilon*(rx**2)))
            Fy += (4*epsilon* ((12*((sigma**12)/(ry**13))) - (6*( (sigma**6) / (ry**7))) )) + ((charge**2) / (varepsilon*(ry**2)))
            Fz += (4*epsilon* ((12*((sigma**12)/(rz**13))) - (6*( (sigma**6) / (rz**7))) )) + ((charge**2) / (varepsilon*(rz**2)))
    return (Fx, Fy, Fz)