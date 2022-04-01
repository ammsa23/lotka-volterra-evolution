
__author__ = "Ammaar A. Saeed"
__email__ = "aasaeed@college.harvard.edu"

'''
Description:
------------
This file contains a command line script that completes 
runs of the gLV model under a stochastic framework as
specified in other parts of the repo. 
The script is versatile: the user can specify a 
variety of the parameters of the simulation, 
including the total number of species S present in the
ecosystem and the method of generating the interaction 
matrix. 
'''

# import necessary modules 
import argparse, os
from cmath import exp

import numpy as np 
import matplotlib.pyplot as plt

from parameters_gLV import *
from plot_gLV import *
from simulate_gLV import * 

def get_args(): 
    '''
    Get the arguments from the command line and parse 
    for the simulation. 

    Parameters: 
    -----------
    None

    Returns:
    --------
    argparse.ArgumentParser - ArgumentParser object storing
    all argument information
    '''

    parser = argparse.ArgumentParser(description="A script\
    that simulates the generalized Lotka-Volterra model\
    under a stochastic framework")
    parser.add_argument("--species", type=int, required=True, 
        help="The total number of species to consider in the\
        ecological network")
    parser.add_argument("--seed", type=int, required=True, 
        help="Random seed integer to for reproducibility of\
        random results")
    parser.add_argument("--population_size", type=int, 
        default=10, help="The initial population size of the\
        interacting species in the ecological network, and\
        total size of species of interest population with one\
        mutant individual")
    parser.add_argument("--nsims", type=int, default=10000, 
        help="The number of simulations of gLV model simulations\
        to complete")
    parser.add_argument("--method", type=str, default="May", 
        help="The parameter generating method for the\
        interaction matrix of the ecological network")
    parser.add_argument("--no_plot", action="store_true",
        help="If true, plotting the results and saving the\
        figures is *skipped!*")
    parser.add_argument("--run_name", type=str,
        default="gLV_run", help="A descriptive prefix for\
        output files")
    parser.add_argument("--out_dir", type=str, 
        help="A directory for storing the simulation results\
        and plots")

    return parser.parse_args()

def generate_model_parameters(
        S: int, 
        seed: int, 
        method = "May"
    ): 
    '''
    Generate the model parameters for the simulation process
    These parameters are generated based on the supplied seed,
    so the values returned will be the same for the same S

    Parameters:
    -----------
    int S - total number of species
    int seed - random seed integer
    string method - the method to generate the interaction matrix:
    the two options are either (1) "May" or (2) "Alessina-Tang"

    Returns:
    --------
    np.array (1) - array of birth and death rates (shape: (S, 2))
    np.array (2) - interaction matrix describing the relationship between 
    the S species (shape: (S+1, S+1))
    '''

    # generate the birth and death rates for each species
    birth_death_rates = generate_birth_death_rates(S = S, seed = seed)
    mutant_birth_death = mutant_birth_death_rate(birth_death = birth_death_rates[-1, :], seed = seed)
    birth_death_rates = np.vstack([birth_death_rates, mutant_birth_death])

    # generate the interaction matrix for each species
    interaction_matrix = generate_interaction_matrix(
        S = S, 
        method = method, 
        seed = seed
    )

    return(
        birth_death_rates, 
        interaction_matrix
    )

def make_interaction_matrix_plot(
        interaction_matrix: float, 
        plot_name: str
    ): 
    '''
    Creates a plot of the interaction matrix and saves to the desired
    directory

    Paraneters:
    -----------
    np.array interaction_matrix - interaction matrix describing the 
    the relationships between all species and the mutant species of
    interest 
    string plot_name - filename for saving plot

    Returns: 
    --------
    None
    '''

    plot_interaction_matrix_CLI(
        interaction_matrix = interaction_matrix, 
        plot_name = plot_name
    )

    return 

def make_birth_death_rates_plot(
        birth_rates: np.array, 
        death_rates: np.array, 
        plot_name: str
    ):
    '''
    Creates a plot of the birth and death rates and saves to the 
    desired directory 

    Parameters:
    -----------
    np.array birth_rates - birth rates of each of the species
    np.array death_rates - death rates of each of the species
    string plot_name - filename for saving plot

    Returns:
    --------
    None
    '''

    plot_birth_death_rates_CLI(
        birth_rates = birth_rates, 
        death_rates = death_rates, 
        plot_name = plot_name
    )

    return 

def simulate_gLV_model(
        nsims: int, 
        S: int, 
        pop_size: int, 
        birth_rates: np.array, 
        death_rates: np.array, 
        interaction_matrix: np.array
    ): 
    '''
    Perform the gLV model simulations using the parameters provided

    Parameters: 
    -----------
    int nsims - total number of simulations
    int S - total number of species
    int pop_size - initial population size of all species in the
    ecological network
    np.array birth_rates - birth rates of each of the species
    np.array death_rates - death rates of each of the species
    np.array interaction_matrix - interaction matrix describing the 
    the relationships between all species and the mutant species of
    interest 

    Returns: 
    --------
    int success_counter - count of number of simulations for which 
    the mutant species of interest fixes
    int reduction_counter - count of number of simulations in which 
    the conditions mirror that of the simple Moran process (i.e., all
    other species except species of interest die out)
    list simulated_lengths - array of simulation length to fixation or 
    extinction
    list exp_wait_times - array of waiting times for each transition
    list wt_states - array of population numbers for each simulation at 
    each transition
    list mutant_states - array of population numbers for each simulation
    at each transition
    '''

    # initialize storage for simulation results 
    success_counter = 0 # counts of mutant fixation 
    reduction_counter = 0 # counts of reduction into classical Moran process
    simulated_lengths = list() # number of transitions to absorption 
    exp_wait_times = list() # length of time to absorption 
    wt_states = list() # states of the wt until absorption 
    mutant_states = list() # states of the mutant until absorption 

    # perform nsims simulations
    for i in range(nsims): 
        simulation_results = run_simulation(
            S = S, 
            pop_size = pop_size, 
            birth_rates = birth_rates, 
            death_rates = death_rates, 
            interaction_matrix = interaction_matrix
        )

        # store the results 
        success_counter += simulation_results[0]
        reduction_counter += np.any(simulation_results[1][:,0] == 0)
        simulated_lengths.append(simulation_results[1].shape[0])
        exp_wait_times.append(simulation_results[2].sum())
        wt_states.append(list(simulation_results[1][:,-2]))
        mutant_states.append(list(simulation_results[1][:,-1]))

    return(
        success_counter, 
        reduction_counter, 
        simulated_lengths, 
        exp_wait_times, 
        wt_states, 
        mutant_states
    )

def report_simulation_statistics(
        success_counter: int, 
        reduction_counter: int, 
        simulated_lengths: list, 
        exp_wait_times: list
    ): 
    '''
    Report quick statistics from the simulation run 

    Parameters: 
    -----------
    int success_counter - count of number of simulations for which 
    the mutant species of interest fixes
    int reduction_counter - count of number of simulations in which 
    the conditions mirror that of the simple Moran process (i.e., all
    other species except species of interest die out)
    list simulated_lengths - array of simulation length to fixation or 
    extinction
    list exp_wait_times - array of waiting times for each transition

    Returns:
    --------
    string - summary statistics of simulation 
    '''

    summary_lines = f"""
The estimated fixation probability under the gLV model is {np.round(success_counter / 10000, 4)}
The gLV model was randomly reduced to the classical Moran process {reduction_counter} times
The average number of transitions until absorption for the species of interst is {np.round(np.mean(simulated_lengths), 4)}
The average total wait time until absorption for the species of interest is {np.round(np.mean(exp_wait_times), 4)}
"""

    print(summary_lines)

    return summary_lines 

def make_mutant_population_trajectories(
        mutant_states: np.array, 
        simulated_lengths: np.array, 
        plot_name: str
    ): 
    '''
    Creates a plot of the mutant population over time for all simulations
    and saves to the desired directory 

    Parameters:
    -----------
    np.array mutant_states - array of mixed length describing the 
    population sizes of the mutant species of interest over time 
    np.array simulated_lengths - array of lengths for each simulation 
    str plot_name - filename for saving plot 

    Returns: 
    --------
    None
    '''

    plot_mutant_trajectories_CLI(
        mutant_states = mutant_states, 
        simulated_lengths = simulated_lengths, 
        plot_name = plot_name
    )

def make_wt_population_trajectories(
        wt_states: np.array, 
        simulated_lengths: np.array, 
        plot_name: str
    ): 
    '''
    Creates a plot of the wild-type populations over time for all 
    simulations and saves to the desired trajectory

    Parameters: 
    -----------
    np.array wt_states - array of mixed length describing the 
    population sizes of the wt species of interest over time 
    np.array simulated_lengths - array of lengths for each simulation 
    str plot_name - filename for saving plot 

    Returns: 
    --------
    None
    '''

    plot_wt_trajectories_CLI(
        wt_states = wt_states, 
        simulated_lengths = simulated_lengths, 
        plot_name = plot_name
    )

    return 

if __name__ == "__main__": 
    '''
    Run the simulation.
    '''

    # parse the arguments 
    args = get_args()

    # generate the model parameters 
    birth_death_rates, interaction_matrix = generate_model_parameters(
        S = args.species, 
        seed = args.seed, 
        method = args.method
    )

    # make the interaction matrix plot
    make_interaction_matrix_plot(
        interaction_matrix = interaction_matrix, 
        plot_name = f"{args.out_dir}/{args.run_name}_interaction_matrix.png"
    )

    # make birth and death rates plot
    make_birth_death_rates_plot(
        birth_rates = birth_death_rates[:, 0], 
        death_rates = birth_death_rates[:, 1], 
        plot_name = f"{args.out_dir}/{args.run_name}_birth_death_rates.png"
    )

    # open a file to store a summary of the simulation results
    f = open(f"{args.out_dir}/{args.run_name}_summary.txt")

    # calculate and report the fixation probability for a simple Moran process
    fixation_probability = calculate_moran_fixation_probability(
        S = args.species, 
        pop_size = args.pop_size, 
        birth_rates = birth_death_rates[:, 0], 
        death_rates = birth_death_rates[:, 1]
    )
    fixation_line = f"The fixation probability of the standard Moran process is \
        {np.round(fixation_probability, 4)}"
    f.write(fixation_line+"\n\n")
    print(fixation_line, end="\n\n")

    # run the gLV model simulation 
    success_counter, \
    reduction_counter, \
    simulated_lengths, \
    exp_wait_times, \
    wt_states, \
    mutant_states = simulate_gLV_model(
        nsims = args.nsims, 
        S = args.species, 
        pop_size = args.pop_size, 
        birth_rates = birth_death_rates[:, 0], 
        death_rates = birth_death_rates[:, 1], 
        interaction_matrix = interaction_matrix
    )

    # report statistics from the simulation and store in summary file
    summary_lines = report_simulation_statistics(
        success_counter = success_counter, 
        reduction_counter = reduction_counter, 
        simulated_lengths = simulated_lengths, 
        exp_wait_times = exp_wait_times
    )
    f.write(summary_lines)
    f.close()

    # make trajectory plots for wt and mutant of species of interest
    make_wt_population_trajectories(
        wt_states = wt_states, 
        simulated_lengths = simulated_lengths, 
        plot_name = f"{args.out_dir}/{args.run_name}_wt_trajectories.png"
    )
    make_mutant_population_trajectories(
        mutant_states = mutant_states, 
        simulated_lengths = simulated_lengths, 
        plot_name = f"{args.out_dir}/{args.run_name}_mutant_trajectories.png"
    )
