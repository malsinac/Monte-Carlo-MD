"""
An object that returns some statistics and analysis
"""
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

class Statistical_analysys(object):
    def __init__(self, States_Library, simulation_name):
        self.States_Library = States_Library
        self.simulation_name = simulation_name

    def RMSD(self):
        _Nparticles = len(self.States_Library[0])
        reference = self.States_Library[0]
        a = 0
        RMSD = [0]
        timesteps = [i for i in range(len(self.States_Library))]
        for i in range(1, len(self.States_Library)):
            for j in range(_Nparticles):
                for k in (0, 1, 2):
                    #print(self.States_Library[i][str(j)][k], reference[str(j)][k])
                    a += (self.States_Library[i][str(j)][k] - reference[str(j)][k])**2
            RMSD.append(np.sqrt(a / _Nparticles))
        fig = plt.figure(figsize=(12,6))
        plt.plot(timesteps, RMSD)
        plt.ylabel('RMSD [A]')
        plt.xlabel('Timestep')
        plt.tight_layout()
        plt.savefig(f'RMSD_{self.simulation_name}.png')




