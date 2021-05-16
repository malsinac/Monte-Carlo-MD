# A Monte Carlo simulator for noble gases
This is a Force-Biased Monte Carlo molecular simulator, suited only for noble gases (in this case, Ar). All the assumptions will be made arround that all particles don't interact between them (i.e no chemical reaction). The basic procedure will be the following:
1.  Generate a random microstate of the sistem.
2.  Calculate the potential energy V(r)*.
3.  Accept the movement if the probability of that microstate is greater than a random number
4.  Repeat on step 1
The way how different parts are implemented will be exposed in the further sections.

## Generating random states
We supose that all particles are ubicated in a box with dimensions Xmax, Ymax, Zmax from the center of it. We also asume that the vectorial space is Cartesian. With this in mind, the new states are generate acording with what it is said in page 417 in [1].

## Calculating the potential energy. Force Field
Because there are not chemical bonds, we only account for non interacting terms. These are the Lennard-Jones terms, with a 12-6 equation, and a Columb electrostatic energetic term. Further work will be done to improve the calculation of energetic terms, since the description of the systems depends directly on how the ForceField is formulated.

## Accepting transition
The Metropolis algorithm is implemented here to accept the transistion state of the system. 

# Future work
The project will follow the next objectives:
- Improve the ForceField.
- Calculate thermodynamic properties.
- Improve the random number generator.
- Allow modification of parameters with a configuration file.
- Use argument parser.
- Improve the random states generator.
- Add some statistics.
- Visualitzation (due creating a web server to display praphics, or making output files suitable for visualitzation programs).
- Allow chemical bond.

# Bibliography
[1]: Andrew R. Leach. Molecular Modelling: Principles and Applications. Second. Essex: Pearson Education, 2001. isbn: 0-582-38210-6.