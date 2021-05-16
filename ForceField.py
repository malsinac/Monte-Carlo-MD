"""
A ForceField for only Ar molecules. Due no chemical bonding, only non-bonding terms are working out. Inspired in the AMBER ForceField, it will use the Lennard-Jones 12-6 for Van der Waals interactions, and columbic interactions.
Ar-Ar parameters -> https://courses.physics.illinois.edu/phys466/sp2013/projects/2001/team5/node13.html
"""
import numpy as np

def PotentialEnergy(positions):
    _Nparticles = len(positions)
    V = 0
    for i in range(0, _Nparticles-1):
        for j in range(i+1, _Nparticles):
            r = np.sqrt((positions[str(i)][0]-positions[str(j)][0])**2 + (positions[str(i)][1]-positions[str(j)][1])**2 + (positions[str(i)][2]-positions[str(j)][2])**2)
            epsilon = 119.8*1.380649*(10**-23) # Ar-Ar value
            sigma = 3.41 #Ar-Ar value
            charge = 1.0 #Atom charge. Must be changed because it is an aproximation
            varepsilon = 1 * 8.8542*10**-12 #Permeability of medium. If first multiplication term is 1, it means we're in vacuum
            V += (4*epsilon*( (sigma/r)**12 - (sigma/r)**6)) + ((charge**2) / (varepsilon*r))       
    return V
