
__author__ = "Ammaar A. Saeed"
__email__ = "aasaeed@college.harvard.edu"

'''
Description: 
------------
This script contains functions for simulating the generalized 
Lotka-Volterra Model with evolution for a total of S species 
under a stochastic framework, particularly a Moran process. 
The functions expect parameters for the simulation (birth and 
death rates, interaction matrix) to be already generated. 

Implementation: 
---------------
The current model is based on a S+1-dimensional 
continuous-time Moran process, with each dimension representing the
population of a certain species. 
We use a finite population Moran process model to simulate 
the potential fixation of rng mutants that arise in the
population for a single mutant, and we repeat these 
simulations many times to estimate the fixation probabilities of
the mutant. 
As such, the simulation uses the embedded chain transition 
matrix to determine the next state of the chain and then 
calculates the exponential waiting time for this step 
accordingly. Since there are technically two steps that 
occur for each transition, we simply assign the time for 
the longer step as the waiting time since the other will 
have already happened by then. 
'''

import numpy as np 
from numpy.random import default_rng

# define a standard rng object for all functions 
rng = default_rng()

def calculate_adjusted_birth_rates(
        system_state: np.array, 
        birth_rates: np.array, 
        interaction_matrix: np.array
    ): 
    '''
    Calculate the adjusted birth rates for each of the species 
    given the current state of the system 

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

    return system_state * birth_rates + system_state @ \
        ((interaction_matrix > 0) * interaction_matrix)


def calculate_adjusted_death_rates(
        system_state: np.array, 
        death_rates: np.array, 
        interaction_matrix: np.array
    ): 
    '''
    Calculate the adjusted death rates for each of the species
    given the current state of the system

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

    return system_state * death_rates + system_state @ \
        ((interaction_matrix > 0) * interaction_matrix)

def calculate_step_probability(
        system_state: np.array, 
        adj_birth_rates: np.array,
        adj_death_rates: np.array
    ): 
    '''
    Calculate the probabilities for each of the potential transitions
    given the current state of the system. 
    The possible transitions require that a single individual in 
    one of the species births/replicates another individual of that
    species, and a single individual in another of the species 
    dies. We treat these two events as indepdent, scaling linearly
    with the population sizes in accordance with minimum 
    Exponential waiting times. 

    Parameters: 
    -----------
    np.array system_state - current populations of all species and 
    mutant species of interest in the system 
    np.array adj_birth_rates - birth rates adjusted by contributions
    from the interaction matrix 
    np.array adj_death_rates - death rates adjusted by contributions 
    from the interaction matrix

    Returns:
    --------
    np.array - array of probabilities for potential transitions; 
    probabilities are arranged in the order of the species given as 
    first birth events for each species and then death events 
    (shape: (2S,))
    '''

    # calculate the birth probabilities 
    birth_probs = system_state * adj_birth_rates / \
        np.sum(system_state * adj_birth_rates)

    # calculate the death probabilities 
    death_probs = system_state * adj_death_rates / \
        np.sum(system_state * adj_death_rates)

    return birth_probs, death_probs

def run_step(
        system_state: np.array, 
        birth_rates: np.array, 
        death_rates: np.array, 
        interaction_matrix: np.array
    ): 
    '''
    Run a single simulation step by taking in the current population 
    of the system, calculating the transition probabilities, 
    rngly determining the next transition, and drawing an 
    exponential waiting time

    Parameters: 
    -----------
    np.array system_state - current populations of all species and 
    mutant species of interest in the system 
    np.array birth_rates - birth rates of each of the species 
    np.array death_rates - death rates of each of the species
    np.array interaction_matrix - interaction matrix describing the 
    the relationships between all species and the mutant species of
    interest 

    Returns: 
    --------
    np.array (1) - populations of all species and mutant species of 
    interest in the system after the current step 
    np.array (2) - waiting time for the process to reach this point 
    '''

    # calculate the adjusted birth and death rates 
    adj_birth_rates = calculate_adjusted_birth_rates(
        system_state, 
        birth_rates, 
        interaction_matrix
    )
    adj_death_rates = calculate_adjusted_death_rates(
        system_state, 
        death_rates, 
        interaction_matrix
    )

    # calculate the transition probabilities according to these rates 
    birth_probs, death_probs = calculate_step_probability(
        system_state, 
        adj_birth_rates, 
        adj_death_rates
    )

    # determine the species replicating in this step 
    unif_sample = rng.uniform(low = 0, high = 1)
    birth_idx = np.nonzero(np.cumsum(birth_probs) > unif_sample)[0]

    # renormalize the death probabilities 
    remove_birth = np.ones(system_state.shape)
    remove_birth[birth_idx] = 0
    death_probs = death_probs * remove_birth
    death_probs = death_probs / np.sum(death_probs)

    # determine the species dying in this step 
    unif_sample = rng.uniform(low = 0, high = 1)
    death_idx = np.nonzero(np.cumsum(death_probs) > unif_sample)[0]

    # simulate the exponential waiting time, taking the max of the minimum of
    # exponential rng variables 
    exp_wait_time = np.max(
        rng.exponential(system_state[birth_idx] * adj_birth_rates[birth_idx]), 
        rng.exponential(system_state[death_idx] * adj_death_rates[death_idx])
    )
    
    # update the populations of the system 
    system_state[birth_idx] += 1
    system_state[death_idx] -= 1

    return system_state, exp_wait_time

def run_simulation(
        S: int, 
        birth_rates: np.array, 
        death_rates: np.array, 
        interaction_matrix: np.array
    ): 
    '''
    Run a simulation of the generalized Lotka-Volterra under the 
    stochastic framework, with stopping conditions of either (i) 
    the species of interest goes extinct or (ii) the mutant 
    species of interest

    Parameters:
    -----------
    int S - total number of species 
    np.array birth_rates - birth rates of each of the species 
    np.array death_rates - death rates of each of the species
    np.array interaction_matrix - interaction matrix describing the 
    the relationships between all species and the mutant species of
    interest 

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
    population0 = 10
    system_state = np.hstack([np.repeat(population0, S), 1])

    # create arrays to store the waiting times and states of the system
    system_states = list()
    exp_wait_times = list()

    # set counter to count total number of steps/transitions made 
    counter = 0

    # use a while loop to determine if a stopping condition is met 
    while (
        system_state[-2] > 0 and 
        system_state[-1] > 0 and 
        counter < 10^5
    ): 
        
        # simulate a single step in the loop 
        system_state, exp_wait_time = run_step(system_state, birth_rates, death_rates, interaction_matrix)

        # append the current information about the system trajectory
        system_states.append(system_state)
        exp_wait_times.append(exp_wait_time)

    # return the collected data for the entire simulation 
    return np.array(system_states), np.array(exp_wait_times)