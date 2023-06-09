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

def test(params):
    D, E, phi = runSimulation('Exercice6_2023_student', params)
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
    plt.plot(phi[:, 0], phi[:, 1], label='Computed $\phi$')
    plt.plot(phi[:, 0], phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va']), label='Analytical $\phi$')
    # plt.rcParams.update({
    # "text.usetex": True,
    # "font.family": "Helvetica"
    # })
    plt.legend()
    plt.xlabel('$r$')
    plt.ylabel('$\phi$')
    plt.show()



# convergence study of phi as a function of N1, compute the error using the analytical solution
def convergence_study_phi(params):
    nb_points = 100
    # nb_points = input('Number of points for convergence study: ')
    N1 = np.logspace(1, 3, nb_points)
    params['N2'] = 1
    error = np.zeros(len(N1))
    for i in range(nb_points):
        params['N1'] = int(N1[i])
        D, E, phi = runSimulation('Exercice6_2023_student', params)
        error[i] = np.max(np.abs(phi[:, 1] - phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va'])))
    # plot error with a log scale on the y axis
    plt.semilogx(N1, error)
    plt.xlabel('N1')
    plt.ylabel('error')
    plt.title('Convergence study of phi')
    plt.show()

def convergence_study_phi_half(params):
    nb_points = 100
    # nb_points = input('Number of points for convergence study: ')
    N = np.logspace(1, 4, nb_points)
    error = np.zeros(len(N))
    for i in range(nb_points):
        params['N1'] = int(N[i])/2
        params['N2'] = int(N[i])/2
        D, E, phi = runSimulation('Exercice6_2023_student', params)
        error[i] = np.max(np.abs(phi[:, 1] - phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va'])))
    plt.loglog(N, error)
    plt.xlabel('N')
    plt.ylabel('error')
    plt.title('Convergence study of phi')
    plt.show()

def convergence_study_phi_final_version(params):
    nb_points = 100
    # nb_points = input('Number of points for convergence study: ')
    N = np.logspace(1, 4, nb_points)

    error1 = np.zeros(len(N))
    params['p'] = 1.0
    for i in range(nb_points):
        params['N1'] = int(N[i])/2
        params['N2'] = int(N[i])/2
        D, E, phi = runSimulation('Exercice6_2023_student', params)
        error1[i] = np.max(np.abs(phi[:, 1] - phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va'])))

    error2 = np.zeros(len(N))
    params['p'] = 0.0
    for i in range(nb_points):
        params['N1'] = int(N[i])/2
        params['N2'] = int(N[i])/2
        D, E, phi = runSimulation('Exercice6_2023_student', params)
        error2[i] = np.max(np.abs(phi[:, 1] - phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va'])))

    error3 = np.zeros(len(N))
    params['p'] = 0.5
    for i in range(nb_points):
        params['N1'] = int(N[i])/2
        params['N2'] = int(N[i])/2
        D, E, phi = runSimulation('Exercice6_2023_student', params)
        error3[i] = np.max(np.abs(phi[:, 1] - phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va'])))

    # plt.rcParams.update({
    #     "text.usetex": True,
    #     "font.family": "Helvetica"
    # })

    # plot errors in loglog scale with different linewidths
    plt.loglog(N, error3, label='$p=0.5$', linewidth=9)
    plt.loglog(N, error1, label='$p=1$', linewidth=3)
    plt.loglog(N, error2, label='$p=0$', linewidth=1)
    plt.xlabel('N')
    plt.ylabel('error')
    # plt.title('Convergence study of $\phi$')
    plt.legend()
    plt.show()



def convergence_study_phi_test(params):
    nb_points = 70
    # nb_points = input('Number of points for convergence study: ')
    N1 = np.logspace(1, 3, nb_points)
    params['N2'] = 1
    error1 = np.zeros(len(N1))
    for i in range(nb_points):
        params['N1'] = int(N1[i])
        D, E, phi = runSimulation('Exercice6_2023_student', params)
        error1[i] = np.max(np.abs(phi[:, 1] - phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va'])))

    params['N2'] = 2
    error2 = np.zeros(len(N1))
    for i in range(nb_points):
        params['N1'] = int(N1[i])
        D, E, phi = runSimulation('Exercice6_2023_student', params)
        error2[i] = np.max(np.abs(phi[:, 1] - phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va'])))

    params['N2'] = 3
    error3 = np.zeros(len(N1))
    for i in range(nb_points):
        params['N1'] = int(N1[i])
        D, E, phi = runSimulation('Exercice6_2023_student', params)
        error3[i] = np.max(np.abs(phi[:, 1] - phi_analytical(phi[:, 0], params['ra'], params['R'], params['Va'])))

    # plot error with a log scale on the y axis
    plt.semilogx(N1, error1)
    plt.semilogx(N1, error2)
    plt.semilogx(N1, error3)
    plt.xlabel('N1')
    plt.ylabel('error')
    plt.title('Convergence study of phi')
    plt.show()






# test()
# basic_plot()

# c)
parameters['ra'] = 0.1
parameters['rb'] = 0.55
parameters['R'] = 1.0
parameters['N1'] = 500
parameters['N2'] = 500
parameters['p'] = 0.5
parameters['talk'] = 'false'
# plot_phi(parameters)
parameters['verbose'] = 0
# convergence_study_phi_half(parameters)
# convergence_study_phi_test(parameters)
# parameters['p'] = 1.0
# plot_phi(parameters)
# convergence_study_phi(parameters)
# parameters['p'] = 0.0
# plot_phi(parameters)
# convergence_study_phi(parame  ters)
# parameters['p'] = 0.5
# convergence_study_phi(parameters)
# convergence_study_phi_final_version(parameters)


# test dimitri

parameters_dimitri = {
    # Physique:
    'ra':0.03,
    'rb':0.05,
    'R':0.10,
    'A':2e4,
    'epsilon_a':1,
    'epsilon_b':5,
    'epsilon_R':2,
    'Va':1.5,
    'VR':0.0,

    # Numérique:
    'N1':10,
    'N2':10,
    'output':'output',
    'p':0.5,

    # Set to 0 to reduce printouts
    'verbose':1 ,
    'talk':'true',
}

test(parameters_dimitri)