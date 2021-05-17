import time
iniit = time.time()
import tracemalloc
tracemalloc.start()
#### Start of the script ####
import numpy as np
import random
from ForceField import PotentialEnergy, Force
import uuid
import datetime

#Phisic Variables
kb = 1.380649*(10**-23) # Boltzman constant
A = 0.5 #For the biased-Monte Carlo force desplacement. Need to confirm that parameter 

#User variables
_Ntimesteps = 500
_Nparticles = 10
_Temperature = 300. #In Kelvins
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
_NStatesAccepted = 0
_NStatesNotAccepted = 0
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
            InitialState[str(j)][0] + (((A * Force(ActualState, NewStateCandidate)[0])/(kb * _Temperature))+np.random.normal(0, 2*A)),
            InitialState[str(j)][1] + (((A * Force(ActualState, NewStateCandidate)[1])/(kb * _Temperature))+np.random.normal(0, 2*A)),
            InitialState[str(j)][2] + (((A * Force(ActualState, NewStateCandidate)[2])/(kb * _Temperature))+np.random.normal(0, 2*A))
            )
        else:
            NewState[str(j)] = ActualState[str(j)]
        print(Force(ActualState, NewStateCandidate)[0], Force(ActualState, NewStateCandidate)[1], Force(ActualState, NewStateCandidate)[2])
    #Compute the potential energy of that microstate
    V_new = PotentialEnergy(NewState)
    #Transition probability
    P = min(
        1,
        np.exp(- (V_new[0] - V_actual[0]) / (kb * _Temperature))
    )
    #Transition procedure
    ActualState_backup = ActualState
    if random.random() <= P:
        transition = 'True'
        ActualState = NewState
        _NStatesAccepted += 1
    else:
        transition = 'False'
        _NStatesNotAccepted += 1

    #Output formating
    with open(f'MonteCarlo_simulation_{simulation_name}.info', 'a+') as output_file:
        output_file.write(
            f"""\n
            Time step: {i}    Particle changed: {str(NewStateCandidate)}    Transition: {transition}
            Actual state -> {ActualState_backup[str(NewStateCandidate)]}
            New state ->    {NewState[str(NewStateCandidate)]}
            Probability of transition: {P}
            Potential energy of actual state -> L-J: {V_actual[1][0]} | Columb: {V_actual[1][1]} | Sum: {V_actual[0]}
            Potential energy of new state ->    L-J: {V_new[1][0]} | Columb: {V_new[1][1]} | Sum: {V_new[0]}
            """
        )
    
#### End of the script ####
maxx, _ = tracemalloc.get_traced_memory()
tracemalloc.stop()
with open(f'MonteCarlo_simulation_{simulation_name}.info', 'a+') as output_file:
    output_file.write(
        f"""
        Simulation done with a total of {_NStatesAccepted} accepted states (proportion of {(_NStatesAccepted) / (_NStatesAccepted+_NStatesNotAccepted)}) and a total of {_NStatesNotAccepted}
        not accepted states (proportion of {(_NStatesNotAccepted) / (_NStatesAccepted+_NStatesNotAccepted)})
        """
    )
    output_file.write(
        """\n
        Total execution time: {:.2f} seconds. Maximum memory allocation: {:.2f} MB
        Seconds/timestep: {:.2f}
        """.format(time.time()-iniit, maxx/(10**6), (time.time()-iniit)/(_Ntimesteps))
    )
print("Total execution time: {:.2f} seconds. Maximum memory allocation: {:.2f} MB".format(time.time()-iniit, maxx/(10**6)))