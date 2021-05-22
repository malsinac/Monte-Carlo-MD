import time
iniit = time.time()
import tracemalloc
tracemalloc.start()
#### Start of the script ####
import numpy as np
import random
from ForceField import PotentialEnergy
from RandomStates import CreateStates
from Statistics import Statistical_analysys
import uuid
import datetime

#Phisic Variables
kb = 0.008314463 # Boltzman constant in kJ/(mol·K)

#User variables
_Ntimesteps = 10000
_Nparticles = 10
_Temperature = 3000 #In Kelvins
_RandomInitial = True

#Box dimensions
_Xmax, _Ymax, _Zmax = 50, 50, 50

#Starting of simulation
simulation_name = uuid.uuid4()
print(f'Simulation name -> {simulation_name}')
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
        InitialState[str(i)] = ((2*random.random()-1)*_Xmax, (2*random.random()-1)*_Ymax, (2*random.random()-1)*_Zmax)
#Output
with open(f'MonteCarlo_simulation_{simulation_name}.info', 'a+') as output_file:
    output_file.write('\nThe initial position of atoms is the following:')
    for key in InitialState.keys():
        output_file.write(f'\nAtom {key} -> {InitialState[key]}')

#Starting the main procedure
_NStatesAccepted = 0
_NStatesNotAccepted = 0
States_Library = [InitialState]
for i in range(_Ntimesteps):
    if i == 0:
        ActualState = InitialState   
    #Calculate the potential energy of actual state
    V_actual = PotentialEnergy(ActualState)
    #Generate a random state
    config_file = {
        'method': 3,
        'A': 0.05,
        'max_displacement': 10
    }
    NewOne = CreateStates(ActualState, config_file, kb=kb, Temp=_Temperature, box_size=(_Xmax, _Ymax, _Zmax), N_acepted=(_NStatesAccepted/(_NStatesAccepted+_NStatesNotAccepted+np.spacing(0))))
    NewState = NewOne[0]
    NewStateCandidate = NewOne[1]
    #Compute the potential energy of that microstate
    V_new = PotentialEnergy(NewState)
    #Transition probability
    P = min(
        1,
        np.exp(- ((V_new[0] - V_actual[0]) / (kb * _Temperature)))
    )
    #Transition procedure
    ActualState_backup = ActualState
    if random.random() < P:
        transition = 'True'
        ActualState = NewState
        _NStatesAccepted += 1
    else:
        transition = 'False'
        _NStatesNotAccepted += 1
    
    States_Library.append(ActualState)
    #Output formating
    with open(f'MonteCarlo_simulation_{simulation_name}.info', 'a+') as output_file:
        output_file.write(
            f"""\n
            Time step: {i}    Particle changed: {str(NewStateCandidate)}    Transition: {transition}
            Actual state -> {ActualState_backup[str(NewStateCandidate)]}
            New state ->    {NewState[str(NewStateCandidate)]}
            Probability of transition: {P}
            Potential energy of actual state -> L-J: {V_actual[1][0]} kJ/mol | Columb: {V_actual[1][1]} kJ/mol | Sum: {V_actual[0]} kJ/mol
            Potential energy of new state ->    L-J: {V_new[1][0]} kJ/mol | Columb: {V_new[1][1]} kJ/mol | Sum: {V_new[0]} kJ/mol
            """
        )

#Analysis of simulation
simu_stats = Statistical_analysys(States_Library, simulation_name)
simu_stats.RMSD()

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