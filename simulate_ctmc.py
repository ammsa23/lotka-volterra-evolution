
__author__ = "Ammaar A. Saeed"
__email__ = "aasaeed@college.harvard.edu"

'''
Description: 
------------
This script contains functions for simulating the generalized 
Lotka-Volterra Model with evolution for a total of S species 
under a stochastic framework. 
The functions expect parameters for the simulation (birth and 
death rates, interaction matrix) to be already generated. 

Implementation: 
---------------
The current model is based on a S+1-dimensional 
continuous-time Markov chain, with each dimension representing the
population of a certain species. 
As such, the simulation uses the embedded chain transition 
matrix to determine the next state of the chain and then 
calculates the exponential waiting time for this step 
accordingly. 
'''

import numpy as np 

def calculate_adjusted_birth_rates(birth_rates: np.array, 
                                   interaction_matrix: np.array): 
    '''
    Calculate the adjusted birth rates for each of the species

    Parameters:
    -----------
    np.array birth_rates - birth rates of each of the species 
    np.array interaction_matrix - interaction matrix describing the 
    the relationships between all species and the mutant species of
    interest 

    Returns: 
    --------
    np.array - array of adjusted birth rates for each of the species, 
    adding only the positive interactions 
    '''

    return birth_rates + \
        np.sum((interaction_matrix > 1) * interaction_matrix, axis = 1)

def calculate_adjusted_death_rates(death_rates: np.array, 
                                   interaction_matrix: np.array): 
    '''
    Calculate the adjusted death rates for each of the species

    Parameters:
    -----------
    np.array death_rates - death rates of each of the species 
    np.array interaction_matrix - interaction matrix describing the 
    the relationships between all species and the mutant species of
    interest 

    Returns: 
    --------
    np.array - array of adjusted death rates for each of the species, 
    adding only the negative interactions 
    '''

    return death_rates - \
        np.sum((interaction_matrix < 1) * interaction_matrix, axis = 1)

def calculate_step_probability(): 
    pass

def run_step(): 
    pass

def run_simulation(S: int, adj_birth_rates: np.array, 
                   adj_death_rates: np.array): 
    '''
    Run a simulation of the generalized Lotka-Volterra under the 
    stochastic framework, with stopping conditions of either (i) 
    the species of interest goes extinct or (ii) the mutant 
    species of interest

    Parameters:
    -----------
    int S - total number of species 
    np.array adj_birth_rates - birth rates adjusted by contributions
    from the interaction matrix 
    np.array adj_death_rates - death rates adjusted by contributions 
    from the interaction matrix

    Returns: 
    --------
    bool - boolean data type that returns whether or not the mutant 
    species of interest fixes (1 if the mutant fixes; 0 if the mutant
    goes extinct)
    np.array (1) - array of exponential wait times between transitions
    (shape: (steps,))
    np.array (2) - array of states describing the population sizes of 
    each species at each step (shape: (steps, 3))
    '''

    # initialize the populations of each of the species 
    populations = np.hstack([np.repeat(10, S), 1])

    # create arrays to store the waiting times and states of the system
    exp_wait_times = list()
    system_states = list()

    # set counter to count total number of steps/transitions made 
    counter = 0

    # use a conditional while loop 
    while (populations[-2] > 0 and populations[-1] > 0 and 
           counter < 10^5): 
        pass 
