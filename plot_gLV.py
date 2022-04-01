
__author__ = "Ammaar A. Saeed"
__email__ = "aasaeed@college.harvard.edu"

'''
Description: 
------------
This file contains a variety of functions meant for plotting 
information and results of the generalized Lotka-Volterra 
Model with evolution simulation studies. 
'''

import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.lines as lines

# set the font parameters 
font = {
    "family": "serif"
    }
plt.rc("font", **font)

#####################################################
############# JUPYTER NOTEBOOK PLOTTING #############
#####################################################

def plot_interaction_matrix(
        interaction_matrix: np.array, 
        cmap = "bwr"
    ): 
    '''
    Plot the interaction matrix for all species and the mutant
    species of interest, labeling the matrix as such
    The mutant is labeled with an asterisk (*)

    Parameters: 
    -----------
    np.array interaction_matrix - interaction matrix describing the 
    the relationships between all species and the mutant species of
    interest 
    str cmap - matplotlib cmap to be used for plotting; default is
    bwr

    Returns: 
    --------
    None
    '''

    # initialize figure and axis objects for the plot 
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))

    # plot the interaction matrix 
    mat = ax.matshow(
        interaction_matrix, 
        cmap = cmap, 
        vmin = -1.25, 
        vmax = 1.25
    )

    # create the tick labels 
    tick_labels = np.array(
        [
            f"S{i+1}" if i != interaction_matrix.shape[0] - 1 else f"S{i}*" \
            for i in np.arange(interaction_matrix.shape[0])
        ]
    )

    # labeling the plot 
    ax.set_title("Species Interaction Matrix", fontsize=20)
    ax.set_xlabel("Species", fontsize=20)
    ax.set_ylabel("Species", fontsize=20)
    ax.set_xticks(np.arange(interaction_matrix.shape[0]))
    ax.set_yticks(np.arange(interaction_matrix.shape[0]))
    ax.set_xticklabels(tick_labels, fontsize=20, rotation="vertical")
    ax.set_yticklabels(tick_labels, fontsize=20)
    ax.xaxis.set_ticks_position('bottom')

    # add a colorbar 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.3)
    cbar = plt.colorbar(mat, cax=cax)
    cbar.ax.tick_params(labelsize=14)

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
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))

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
    ax.set_title("Birth/Death Rate of Species", fontsize=20)
    ax.set_xlabel("Species", fontsize=20)
    ax.set_ylabel("Birth/Death Rate", fontsize=20)
    ax.set_xticks(np.arange(birth_rates.shape[0]))
    rate_max = np.max(np.hstack([birth_rates, death_rates]))
    ax.set_yticks(np.arange(0, (rate_max // 0.2 + 1) * 0.2 + 0.1, 0.2))
    ax.set_xticklabels(tick_labels, fontsize=20, rotation="vertical")
    ax.set_yticklabels(np.round(ax.get_yticks(), 2), fontsize=20)

    # make a legend
    ax.legend(fontsize=16)

    # show the figure
    plt.show()

def plot_mutant_trajectories(
        mutant_states: np.array, 
        simulated_lengths: np.array
    ): 
    '''
    Plot the mutant population over time for all simulations 

    Parameters: 
    -----------
    np.array mutant_states - array of mixed length describing the 
    population sizes of the mutant species of interest over time 
    np.array simulated_lengths - array of lengths for each simulation 

    Returns: 
    --------
    None
    '''

    # initialize figure and axis objects for the plot 
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))

    # plot the mutant trajectories over time, using maximum time for x-axis
    max_time = np.max(simulated_lengths)
    for sim in mutant_states: 

        # adjust sim so that it matches the length of max_time 
        sim = sim + (max_time - len(sim)) * [sim[-1]]

        # color the trajectories based on fixation or extinction 
        if sim[-1] == 0: 
            ax.plot(np.arange(1, max_time + 1), sim, "r-", alpha=0.15)
        else: 
            ax.plot(np.arange(1, max_time + 1), sim, "b-", alpha=0.15)

    # labeling the plot 
    ax.set_title("Mutant Population Trajectories", fontsize=20)
    ax.set_xlabel("Transitions", fontsize=20)
    ax.set_ylabel("Mutant Population", fontsize=20)
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())
    ax.set_xticklabels(ax.get_xticks().astype("int32"), fontsize=20)
    ax.set_yticklabels(ax.get_yticks().astype("int32"), fontsize=20)

    # make custom figure handles 
    extinct = lines.Line2D([], [], color="r", label = "Extinction")
    fixed = lines.Line2D([], [], color="b", label = "Fixation")

    # make a legend 
    ax.legend(
        handles = [extinct, fixed], 
        loc = "upper left", 
        fontsize=16
    )

    # show the figure
    plt.show()

def plot_wt_trajectories(
        wt_states: np.array, 
        simulated_lengths: np.array
    ): 
    '''
    Plot the wt population over time for all simulations 

    Parameters: 
    -----------
    np.array wt_states - array of mixed length describing the 
    population sizes of the wt species of interest over time 
    np.array simulated_lengths - array of lengths for each simulation 

    Returns: 
    --------
    None
    '''

    # initialize figure and axis objects for the plot 
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))

    # plot the wt trajectories over time, using maximum time for x-axis
    max_time = np.max(simulated_lengths)
    for sim in wt_states: 

        # adjust sim so that it matches the length of max_time 
        sim = sim + (max_time - len(sim)) * [sim[-1]]

        # color the trajectories based on fixation or extinction 
        if sim[-1] == 0: 
            ax.plot(np.arange(1, max_time + 1), sim, "r-", alpha=0.15)
        else: 
            ax.plot(np.arange(1, max_time + 1), sim, "b-", alpha=0.15)

    # labeling the plot 
    ax.set_title("Wild-type Population Trajectories", fontsize=20)
    ax.set_xlabel("Transitions", fontsize=20)
    ax.set_ylabel("Wild-Type Population", fontsize=20)
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())
    ax.set_xticklabels(ax.get_xticks().astype("int32"), fontsize=20)
    ax.set_yticklabels(ax.get_yticks().astype("int32"), fontsize=20)

    # make custom figure handles 
    extinct = lines.Line2D([], [], color="r", label = "Extinction")
    fixed = lines.Line2D([], [], color="b", label = "Fixation")

    # make a legend 
    ax.legend(
        handles = [extinct, fixed], 
        loc = "upper left", 
        fontsize=16
    )

    # show the figure
    plt.show()

#####################################################
############### COMMAND LINE PLOTTING ###############
#####################################################

def plot_interaction_matrix_CLI(
        interaction_matrix: np.array, 
        plot_name: str,
        cmap = "bwr", 
        no_plot = False
    ): 
    '''
    Plot the interaction matrix for all species and the mutant
    species of interest, labeling the matrix as such
    The mutant is labeled with an asterisk (*)

    Parameters: 
    -----------
    np.array interaction_matrix - interaction matrix describing the 
    the relationships between all species and the mutant species of
    interest 
    str cmap - matplotlib cmap to be used for plotting; default is
    bwr
    str plot_name - filename for saving plot 
    bool no_plot - if True, then no plots are shown

    Returns: 
    --------
    None
    '''

    # initialize figure and axis objects for the plot 
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))

    # plot the interaction matrix 
    mat = ax.matshow(
        interaction_matrix, 
        cmap = cmap, 
        vmin = -1.25, 
        vmax = 1.25
    )

    # create the tick labels 
    tick_labels = np.array(
        [
            f"S{i+1}" if i != interaction_matrix.shape[0] - 1 else f"S{i}*" \
            for i in np.arange(interaction_matrix.shape[0])
        ]
    )

    # labeling the plot 
    ax.set_title("Species Interaction Matrix", fontsize=20)
    ax.set_xlabel("Species", fontsize=20)
    ax.set_ylabel("Species", fontsize=20)
    ax.set_xticks(np.arange(interaction_matrix.shape[0]))
    ax.set_yticks(np.arange(interaction_matrix.shape[0]))
    ax.set_xticklabels(tick_labels, fontsize=20, rotation="vertical")
    ax.set_yticklabels(tick_labels, fontsize=20)
    ax.xaxis.set_ticks_position('bottom')

    # add a colorbar 
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.3)
    cbar = plt.colorbar(mat, cax=cax)
    cbar.ax.tick_params(labelsize=14)

    # show the figure
    fig.savefig(plot_name, dpi = 200)

    if not no_plot: 
        plt.show()

    return 

def plot_birth_death_rates_CLI(
        birth_rates: np.array, 
        death_rates: np.array, 
        plot_name: str, 
        no_plot = False
    ): 
    '''
    Plot the birth and death rates for all species and mutant 
    species of interest

    Parameters: 
    -----------
    np.array birth_rates - birth rates of each of the species
    np.array death_rates - death rates of each of the species 
    string plot_name - filename for saving plot 
    bool no_plot - if True, then no plots are shown

    Returns: 
    --------
    None
    '''

    # initialize figure and axis objects for the plot 
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))

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
    ax.set_title("Birth/Death Rate of Species", fontsize=20)
    ax.set_xlabel("Species", fontsize=20)
    ax.set_ylabel("Birth/Death Rate", fontsize=20)
    ax.set_xticks(np.arange(birth_rates.shape[0]))
    rate_max = np.max(np.hstack([birth_rates, death_rates]))
    ax.set_yticks(np.arange(0, (rate_max // 0.2 + 1) * 0.2 + 0.1, 0.2))
    ax.set_xticklabels(tick_labels, fontsize=20, rotation="vertical")
    ax.set_yticklabels(np.round(ax.get_yticks(), 2), fontsize=20)

    # make a legend
    ax.legend(fontsize=16)

    # show the figure
    fig.savefig(plot_name, dpi = 200)

    if not no_plot: 
        plt.show()

    return 

def plot_mutant_trajectories_CLI(
        mutant_states: np.array, 
        simulated_lengths: np.array, 
        plot_name: str, 
        no_plot = False
    ): 
    '''
    Plot the mutant population over time for all simulations 

    Parameters: 
    -----------
    np.array mutant_states - array of mixed length describing the 
    population sizes of the mutant species of interest over time 
    np.array simulated_lengths - array of lengths for each simulation 
    str plot_name - filename for saving plot 
    bool no_plot - if True, then no plots are shown

    Returns: 
    --------
    None
    '''

    # initialize figure and axis objects for the plot 
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))

    # plot the mutant trajectories over time, using maximum time for x-axis
    max_time = np.max(simulated_lengths)
    for sim in mutant_states: 

        # adjust sim so that it matches the length of max_time 
        sim = sim + (max_time - len(sim)) * [sim[-1]]

        # color the trajectories based on fixation or extinction 
        if sim[-1] == 0: 
            ax.plot(np.arange(1, max_time + 1), sim, "r-", alpha=0.15)
        else: 
            ax.plot(np.arange(1, max_time + 1), sim, "b-", alpha=0.15)

    # labeling the plot 
    ax.set_title("Mutant Population Trajectories", fontsize=20)
    ax.set_xlabel("Transitions", fontsize=20)
    ax.set_ylabel("Mutant Population", fontsize=20)
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())
    ax.set_xticklabels(ax.get_xticks().astype("int32"), fontsize=20)
    ax.set_yticklabels(ax.get_yticks().astype("int32"), fontsize=20)

    # make custom figure handles 
    extinct = lines.Line2D([], [], color="r", label = "Extinction")
    fixed = lines.Line2D([], [], color="b", label = "Fixation")

    # make a legend 
    ax.legend(
        handles = [extinct, fixed], 
        loc = "upper left", 
        fontsize=16
    )

    # show the figure
    fig.savefig(plot_name, dpi = 200)

    if not no_plot:
        plt.show()

    return 

def plot_wt_trajectories_CLI(
        wt_states: np.array, 
        simulated_lengths: np.array, 
        plot_name: str, 
        no_plot = False
    ): 
    '''
    Plot the wt population over time for all simulations 

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

    # initialize figure and axis objects for the plot 
    fig, ax = plt.subplots(nrows = 1, ncols = 1, figsize=(10,10))

    # plot the wt trajectories over time, using maximum time for x-axis
    max_time = np.max(simulated_lengths)
    for sim in wt_states: 

        # adjust sim so that it matches the length of max_time 
        sim = sim + (max_time - len(sim)) * [sim[-1]]

        # color the trajectories based on fixation or extinction 
        if sim[-1] == 0: 
            ax.plot(np.arange(1, max_time + 1), sim, "r-", alpha=0.15)
        else: 
            ax.plot(np.arange(1, max_time + 1), sim, "b-", alpha=0.15)

    # labeling the plot 
    ax.set_title("Wild-type Population Trajectories", fontsize=20)
    ax.set_xlabel("Transitions", fontsize=20)
    ax.set_ylabel("Wild-Type Population", fontsize=20)
    ax.set_xticks(ax.get_xticks())
    ax.set_yticks(ax.get_yticks())
    ax.set_xticklabels(ax.get_xticks().astype("int32"), fontsize=20)
    ax.set_yticklabels(ax.get_yticks().astype("int32"), fontsize=20)

    # make custom figure handles 
    extinct = lines.Line2D([], [], color="r", label = "Extinction")
    fixed = lines.Line2D([], [], color="b", label = "Fixation")

    # make a legend 
    ax.legend(
        handles = [extinct, fixed], 
        loc = "upper left", 
        fontsize=16
    )

    # show the figure
    fig.savefig(plot_name, dpi = 200)

    if not no_plot: 
        plt.show()

    return 