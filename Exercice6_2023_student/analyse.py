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
    'rb':0.099,
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
    'p':0.5,

    # Set to 0 to reduce printouts
    'verbose':1 ,
    'talk':'true',
}

    

#analytical solutions

# def D_analytical(R, ra, rb, A, epsilon_a, epsilon_b, epsilon_R, Va, VR):

# def E_analytical(R, ra, rb, A, epsilon_a, epsilon_b, epsilon_R, Va, VR):

def phi_analytical(r, ra, R, Va):
    return Va*np.log(r/R)/np.log(ra/R)

def phi_analytical_d(r, ra, R, Va):
    return Va*r/(ra - R) + Va*R/(R - ra)

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
        plt.plot(phi[:, 0], phi_analytical(phi[:, 0], parameters['ra'], parameters['R'], parameters['Va']), label='phi_analytical')
        plt.legend('phi', 'phi_analytical')
        plt.xlabel('r')
        plt.title('phi')
        plt.show()

def plot_r(params):
    D, E, phi = runSimulation('Exercice6_2023_student', params)
    plt.plot(D[:, 0], label='r')
    plt.legend('r')
    plt.xlabel('r')

# c) compare the numerical solution with the analytical solution

def plot_phi(params):
    D, E, phi = runSimulation('Exercice6_2023_student', params)
    plt.plot(phi[:, 0], phi[:, 1], label='phi')
    plt.plot(phi[:, 0], phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va']), label='phi_analytical')
    plt.legend('phi', 'phi_analytical')
    plt.xlabel('r')
    plt.title('phi')
    plt.show()



# convergence study of phi as a function of N1, compute the error using the analytical solution
def convergence_study_phi(params):
    N1 = np.logspace(1, 3, 10)
    params['N2'] = 0
    error = np.zeros(len(N1))
    for i in range(len(N1)):
        params['N1'] = int(N1[i])
        phi = runSimulation('Exercice6_2023_student', params)
        error[i] = np.max(np.abs(phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va']) - phi[:, 1]))
    plt.loglog(N1, error)
    plt.xlabel('N1')
    plt.ylabel('error')
    plt.title('Convergence study of phi')
    plt.show()







# test()
# basic_plot()

# c)
parameters['N1'] = 1000
parameters['N2'] = 1
parameters['p'] = 1.0
parameters['talk'] = 'false'
plot_phi(parameters)
# convergence_study_phi(parameters)