
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

    return parser 

def generate_model_parameters(): 
    '''
    
    '''

    pass

def make_interaction_matrix_plot(): 
    '''
    
    '''

    pass

def make_birth_death_rates_plot():
    '''
    
    '''

    pass


def simulate_gLV_model(
        nsims: int, 
        S: int, 
        pop_size: int, 
        birth_rates: np.array, 
        death_rates: np.array, 
        interaction_matrix: np.array
    ): 
    '''
    
    '''

    # initialize storage for simulation results 
    success_counter = 0 # counts of mutant fixation 
    reduction_counter = 0 # counts of reduction into classical Moran process
    sim_lens = list() # number of transitions to absorption 
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
        sim_lens.append(simulation_results[1].shape[0])
        exp_wait_times.append(simulation_results[2].sum())
        wt_states.append(list(simulation_results[1][:,-2]))
        mutant_states.append(list(simulation_results[1][:,-1]))

    return(
        success_counter, 
        reduction_counter, 
        sim_lens, 
        exp_wait_times, 
        wt_states, 
        mutant_states
    )

def report_simulation_statistics(): 
    '''
    
    '''

    pass

def make_mutant_population_trajectories(): 
    '''
    
    '''

    pass

def make_wt_population_trajectories(): 
    '''
    
    '''

    pass

if __name__ == "__main__": 
    '''
    Run the simulation.
    '''

    # parse the arguments 
    parser = get_args()

    # generate the model parameters 
    generate_model_parameters()

    # make the interaction matrix plot
    make_interaction_matrix_plot()

    # make birth and death rates plot
    make_birth_death_rates_plot()

    # calculate and report the fixation probability for a simple Moran process
    calculate_moran_fixation_probability()

    # run the gLV model simulation 
    simulate_gLV_model()

    # report statistics from the simulation 
    report_simulation_statistics()

    # make trajectory plots for wt and mutant of species of interest
    make_mutant_population_trajectories()
    make_wt_population_trajectories()
