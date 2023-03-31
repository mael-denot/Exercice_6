import os
import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm


# Function to run c++ file for a set of parameters
def runSimulation(cppFile, params):
    '''
    Runs c++ simulation on input dictionary params and returns results as a matrix.
    Creates and then deletes two temporary files to pass parameters into and output data out of the c++ code.
    Results matrix contains output at each timestep as each row. Each variable is it's own column.

    Inputs:
        cppFile: String,    Name of c++ file
        params: dict,   Simulation Parameters

    Outputs:
        simulationData: np.array(nsteps+1, k),  matrix containing all output data across k variables
    '''

    # Set name of temporary output data file
    tempDataFileName = 'output_temp'
    params['output'] = tempDataFileName

    # Create temporary file of parameters to get passed into simulation. Deleted after use.
    tempParamsFileName = 'configuration_temp.in'
    with open(tempParamsFileName, 'w') as temp:
        for key in params.keys():
            temp.write("%s %s %s\n" % (key, "=", params[key]))

    # Run simulation from terminal, passing in parameters via temporary file
    os.system(' .\\' + cppFile + ".exe .\\" + tempParamsFileName)

    # Save output data in simData matrix
    simulationData_D   = np.loadtxt(tempDataFileName + '_D.out')
    simulationData_E   = np.loadtxt(tempDataFileName + '_E.out')
    simulationData_phi = np.loadtxt(tempDataFileName + '_phi.out')

    # Delete temp files. Comment out if you want to keep them.
    os.remove(tempParamsFileName)
    os.remove(tempDataFileName + '_D.out')
    os.remove(tempDataFileName + '_E.out')
    os.remove(tempDataFileName + '_phi.out')

    return simulationData_D, simulationData_E, simulationData_phi


parameters = {
    # Physique:
    'ra':0.03,
    'rb':0.05,
    'R':0.10,
    'A':0,
    'epsilon_a':1,
    'epsilon_b':1,
    'epsilon_R':1,
    'Va':1.5,
    'VR':0.0,

    # Numérique:
    'N1':100,
    'N2':100,
    'output':'output',

    # Set to 0 to reduce printouts
    'verbose':1 ,
}

    

#analytical solutions

# def D_analytical(R, ra, rb, A, epsilon_a, epsilon_b, epsilon_R, Va, VR):

# def E_analytical(R, ra, rb, A, epsilon_a, epsilon_b, epsilon_R, Va, VR):

# def phi_analytical(R, ra, rb, A, epsilon_a, epsilon_b, epsilon_R, Va, VR):

# plots

def test():
    D, E, phi = runSimulation('Exercice6_2023_student', parameters)
    print(D, E, phi)

def basic_plot():
    D, E, phi = runSimulation('Exercice6_2023_student', parameters)

    # ask whether to show or not
    show = input('Show plot for D? (y/n) ')
    if show == 'y':
        plt.plot(D[:, 0], D[:, 1], label='D')
        plt.legend('D')
        plt.xlabel('r')
        plt.title('D')
        plt.show()

    show = input('Show plot for E? (y/n) ')
    if show == 'y':
        plt.plot(E[:, 0], E[:, 1], label='E')
        plt.legend('E')
        plt.xlabel('r')
        plt.title('E')
        plt.show()

    show = input('Show plot for phi? (y/n) ')
    if show == 'y':
        plt.plot(phi[:, 0], phi[:, 1], label='phi')
        plt.legend('phi')
        plt.xlabel('r')
        plt.title('phi')
        plt.show()

# c) compare the numerical solution with the analytical solution

# convergence study of phi as a function of N1
def convergence_study_phi():
    N1 = [3, 4, 5, 6, 7, 8, 9, 10]
    phi = []
    for i in range(len(N1)):
        parameters['N1'] = N1[i]
        D, E, phi_i = runSimulation('Exercice6_2023_student', parameters)
        phi.append(phi_i[-1, 1])
    plt.plot(N1, phi, label='phi')
    plt.legend()
    plt.show()

# test()
basic_plot()