"""
A group of functions to build up random micro-states of the system
"""
import random
import numpy as np
from ForceField import Force


def CreateStates(ActualState, config_file, kb, Temp, box_size, N_acepted):
    #It gets a dictionary with the position of actual state, and creates a new one based on it
    _Nparticles = len(ActualState)
    NewState = {}
    NewStateCandidate = random.randint(0, _Nparticles-1) # Choose only one atom to change its position
    for j in range(_Nparticles):
        if j == NewStateCandidate:
            if config_file['method'] == 1:
                NewState[str(j)] = (
                    float(ActualState[str(j)][0] + ((2*random.random() - 1)*config_file['max_displacement'])) % box_size[0],
                    float(ActualState[str(j)][1] + ((2*random.random() - 1)*config_file['max_displacement'])) % box_size[1],
                    float(ActualState[str(j)][2] + ((2*random.random() - 1)*config_file['max_displacement'])) % box_size[2],
                )
            elif config_file['method'] == 2:
                NewState[str(j)] = (
                    float(ActualState[str(j)][0] + (((config_file['A'] * Force(ActualState, NewStateCandidate)[0])/(kb * Temp))+np.random.normal(0, 2*config_file['A']))) % box_size[0],
                    float(ActualState[str(j)][1] + (((config_file['A'] * Force(ActualState, NewStateCandidate)[1])/(kb * Temp))+np.random.normal(0, 2*config_file['A']))) % box_size[1],
                    float(ActualState[str(j)][2] + (((config_file['A'] * Force(ActualState, NewStateCandidate)[2])/(kb * Temp))+np.random.normal(0, 2*config_file['A']))) % box_size[2]
                )
            elif config_file['method'] == 3:
                NewState[str(j)] = (
                    float(ActualState[str(j)][0] + ((2*random.random() - 1)*(np.abs(box_size[0]-ActualState[str(j)][0])*np.exp(1 - (1/(N_acepted**2)))))) % box_size[0],
                    float(ActualState[str(j)][1] + ((2*random.random() - 1)*(np.abs(box_size[1]-ActualState[str(j)][1])*np.exp(1 - (1/(N_acepted**2)))))) % box_size[1],
                    float(ActualState[str(j)][2] + ((2*random.random() - 1)*(np.abs(box_size[2]-ActualState[str(j)][2])*np.exp(1 - (1/(N_acepted**2)))))) % box_size[2],
                )    
        else:
          NewState[str(j)] = ActualState[str(j)]  
    
    return (NewState, NewStateCandidate)