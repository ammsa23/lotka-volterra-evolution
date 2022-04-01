
__author__ = "Ammaar A. Saeed"
__email__ = "aasaeed@college.harvard.edu"

'''
Description: 
------------
This script contains a variety of functions meant for plotting 
information and results of the generalized Lotka-Volterra 
Model with evolution simulation studies. 
'''

import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

# set the font parameters 
font = {
    "family": "monospace"
    }
plt.rc("font", **font)

def plot_interaction_matrix(
        interaction_matrix: np.array, 
        cmap = "cividis"
    ): 
    '''
    Plot the interaction matrix for all species and the mutant
    species of interest, labeling the matrix as such 

    Parameters: 
    -----------
    np.array interaction_matrix - interaction matrix describing the 
    the relationships between all species and the mutant species of
    interest 
    str cmap - matplotlib cmap to be used for plotting; default is
    cividis


    Returns: 
    --------
    None
    '''

    # initialize figure and axis objects for the plot 
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7,7))

    # plot the interaction matrix 
    mat = ax.matshow(
        interaction_matrix, 
        cmap = cmap
    )

    # create the tick labels 
    tick_labels = np.array(
        [
            f"S{i+1}" if i != interaction_matrix.shape[0] - 1 else f"S{i}*" \
            for i in np.arange(interaction_matrix.shape[0])
        ]
    )

    # labeling the plot 
    ax.set_xlabel("Species", fontsize=20)
    ax.set_ylabel("Species", fontsize=20)
    ax.set_xticks(np.arange(interaction_matrix.shape[0]))
    ax.set_yticks(np.arange(interaction_matrix.shape[0]))
    ax.set_xticklabels(tick_labels, fontsize=20)
    ax.set_yticklabels(tick_labels, fontsize=20)
    ax.xaxis.set_ticks_position('bottom')

    # add a colorbar 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.3)
    plt.colorbar(mat, cax=cax)

    # show the figure 
    plt.show()

def plot_birth_death_rates(
        birth_rates: np.array, 
        death_rates: np.array
    ): 
    '''
    Plot the birth and death rates for all species and mutant 
    species of interest

    Parameters: 
    -----------
    np.array birth_rates - birth rates of each of the species
    np.array death_rates - death rates of each of the species 

    Returns: 
    --------
    None
    '''

    # initialize figure and axis objects for the plot 
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(7,7))

    # plot the birth and death rates 
    ax.plot(np.arange(birth_rates.shape[0]), birth_rates, '^-', label = "birth rates")
    ax.plot(np.arange(death_rates.shape[0]), death_rates, 'o-', label = "death rates")

    # create the tick labels 
    tick_labels = np.array(
        [
            f"S{i+1}" if i != birth_rates.shape[0] - 1 else f"S{i}*" \
            for i in np.arange(birth_rates.shape[0])
        ]
    )

    # labeling the plot 
    ax.set_xlabel("Species", fontsize=20)
    ax.set_ylabel("Birth/Death Rate", fontsize=20)
    ax.set_xticks(np.arange(birth_rates.shape[0]))
    ax.set_xticklabels(tick_labels, fontsize=20)
    ax.set_yticklabels(np.round(ax.get_yticks(), 2), fontsize=20)

    # make a legend
    ax.legend()

    # show the figure
    plt.show()

def plot_mutant_trajectories(): 
    '''

    '''

    pass

def plot_trajectory(): 
    '''
    
    '''

    pass