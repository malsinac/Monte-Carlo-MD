import time
iniit = time.time()
import tracemalloc
tracemalloc.start()
#### Start of the script ####
import numpy as np
import random
from ForceField import PotentialEnergy
import uuid
import datetime

#Phisic Variables
kb = 1.380649*(10**-23) # Boltzman constant

#User variables
_Ntimesteps = 100
_Nparticles = 4
_Temperature = 280. #In Kelvins
_RandomInitial = True

#Box dimensions
_Xmax, _Ymax, _Zmax = 10, 10, 10

#Starting of simulation
simulation_name = uuid.uuid4()
with open(f'MonteCarlo_simulation_{simulation_name}.info', 'w') as output_file:
    output_file.write(
        f"""
        Wellcome to MonteCarlo simulation, made with python! This consist on a simulation of {_Nparticles} atoms
        of Argon, in a box of dimensions ({_Xmax}, {_Ymax}, {_Zmax}) Angstroms. The simulation started on {datetime.datetime.now()} 
        """
    )

#Generating an initial random position
if _RandomInitial:
    InitialState = {}
    for i in range(_Nparticles):
        InitialState[str(i)] = ((2*random.random()-1)*_Xmax/2, (2*random.random()-1)*_Ymax/2, (2*random.random()-1)*_Zmax/2)
#Output
with open(f'MonteCarlo_simulation_{simulation_name}.info', 'a+') as output_file:
    output_file.write('\nThe initial position of atoms is the following:')
    for key in InitialState.keys():
        output_file.write(f'\nAtom {key} -> {InitialState[key]}')

#Starting the main procedure
for i in range(_Ntimesteps):
    if i == 0:
        ActualState = InitialState
    #Calculate the potential energy of actual state
    V_actual = PotentialEnergy(ActualState)
    #Generate a random state
    NewState = {}
    NewStateCandidate = random.randint(0, _Nparticles-1) # Choose only one atom to change its position
    for j in range(_Nparticles):
        if j == NewStateCandidate:
            NewState[str(j)] = (
            InitialState[str(j)][0] + (2*random.random()-1)*2,
            InitialState[str(j)][1] + (2*random.random()-1)*2,
            InitialState[str(j)][2] + (2*random.random()-1)*2
        )
        else:
            NewState[str(j)] = ActualState[str(j)]
    #Compute the potential energy of that microstate
    V_new = PotentialEnergy(NewState)
    #Transition probability
    P = min(
        1,
        np.exp(- (V_new - V_actual) / (kb * _Temperature))
    )
    #Transition procedure
    ActualState_backup = ActualState
    if random.random() <= P:
        transition = 'True'
        ActualState = NewState
    else:
        transition = 'False'

    #Output formating
    with open(f'MonteCarlo_simulation_{simulation_name}.info', 'a+') as output_file:
        output_file.write(
            f"""\n
            Time step: {i}    Particle changed: {str(NewStateCandidate)}    Transition: {transition}
            Actual state -> {ActualState_backup[str(NewStateCandidate)]}
            New state ->    {NewState[str(NewStateCandidate)]}
            Probability of transition: {P}
            """
        )
    
#### End of the script ####
maxx, _ = tracemalloc.get_traced_memory()
tracemalloc.stop()
with open(f'MonteCarlo_simulation_{simulation_name}.info', 'a+') as output_file:
    output_file.write("\nTotal execution time: {:.2f} seconds. Maximum memory allocation: {:.2f} MB".format(time.time()-iniit, maxx/(10**6)))
print("Total execution time: {:.2f} seconds. Maximum memory allocation: {:.2f} MB".format(time.time()-iniit, maxx/(10**6)))