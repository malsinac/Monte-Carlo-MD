"""
A ForceField for only Ar molecules. Due no chemical bonding, only non-bonding terms are working out. Inspired in the AMBER ForceField, it will use the Lennard-Jones 12-6 for Van der Waals interactions, and columbic interactions.
Ar-Ar parameters -> https://courses.physics.illinois.edu/phys466/sp2013/projects/2001/team5/node13.html
"""
import numpy as np

epsilon = 0.997 # Ar-Ar value in kJ/mol
sigma = 3.405 #Ar-Ar value in A
charge = -1.358004793 ** -19 #Charge of oxygen in C meanwhile I work out in AR particle charge
varepsilon = 1 / (4 * np.pi * (1.470795318**-30))

#Return potential energy in Joules and force in Newtons

def PotentialEnergy(positions):
    global epsilon, sigma, charge, varepsilon
    _Nparticles = len(positions)
    V_LJ, V_C = 0, 0
    for i in range(0, _Nparticles-1):
        for j in range(i+1, _Nparticles):
            r = np.sqrt((positions[str(i)][0]-positions[str(j)][0])**2 + (positions[str(i)][1]-positions[str(j)][1])**2 + (positions[str(i)][2]-positions[str(j)][2])**2)
            try:
                V_LJ += 4*epsilon*( ((sigma/r)**12) - ((sigma/r)**6))
                V_C  += (varepsilon*(charge**2)) / (r)
            except RuntimeError:
                return 0.
            else:      
                return (V_LJ + V_C, (V_LJ, V_C))

def Force(positions, candidate):
    #Force that actues over a candidate particule
    global epsilon, sigma, charge, varepsilon
    _Nparticles = len(positions)
    Fx_LJ, Fx_C = 0, 0
    Fy_LJ, Fy_C = 0, 0
    Fz_LJ, Fz_C = 0, 0

    for i in range(0, _Nparticles):
        if i == candidate:
            continue
        else:
            rx = np.sqrt((positions[str(i)][0]-positions[str(candidate)][0])**2)
            ry = np.sqrt((positions[str(i)][1]-positions[str(candidate)][1])**2)
            rz = np.sqrt((positions[str(i)][2]-positions[str(candidate)][2])**2)
            try:
                #x-axis
                Fx_LJ += (4*epsilon* ((12*((sigma**12)/(rx**13))) - (6*( (sigma**6) / (rx**7))) ))
                Fx_C  += ((varepsilon*charge**2) / ((rx**2)))
                #y-axis
                Fy_LJ += (4*epsilon* ((12*((sigma**12)/(ry**13))) - (6*( (sigma**6) / (ry**7))) ))
                Fy_C  += ((varepsilon*charge**2) / ((ry**2)))
                #z-axis
                Fz_LJ += (4*epsilon* ((12*((sigma**12)/(rz**13))) - (6*( (sigma**6) / (rz**7))) ))
                Fz_C  += ((varepsilon*charge**2) / ((rz**2)))
            except RuntimeError:
                Fx_LJ, Fx_C, Fy_LJ, Fy_C, Fz_LJ, Fz_C = 0., 0., 0., 0., 0., 0.

    return (Fx_LJ+Fx_C, Fy_LJ+Fy_C, Fz_LJ+Fz_C)