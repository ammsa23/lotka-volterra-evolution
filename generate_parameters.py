
__author__ = "Ammaar A. Saeed"
__email__ = "aasaeed@college.harvard.edu"

'''
Description: 
------------
This script contains functions that generate parameters for the generalized 
Lotka-Volterra Model with evolution represented under a stochastic framework
for S total species and a mutant of one of the species

Parameter Generation: 
---------------------
The birth and death rate parameters are drawn from a folded standard Normal 
distribution 
The interaction matrix parameters are drawn from a standard Normal distribution
The perturbations on the mutant are drawn from a Uniform distribution on 
the interval [-0.1,0.1]
'''

import numpy as np
import numpy.random as random

def generate_birth_death_rates(S: int): 
    '''
    Generate a pair of positive birth and death rates for each species
    
    Parameters: 
    -----------
    int S - total number of species; the number of birth and death rates to 
    generate

    Returns: 
    --------
    np.array - array of birth (column 1) and death (column 2) rates for each 
    species (shape: (S, 2))
    '''

    return np.abs(random.normal(loc = 0, scale = 1, size = (S, 2)))

def mutant_birth_death_rate(birth_death: np.array): 
    '''
    Perturb the given birth and death rates from draws of a Uniform distribution 
    on the interval [-0.1, 0.1]

    Parameters: 
    -----------
    np.array birth_death - birth and death rates of the species of interest; 
    birth and death rates to be perturbed 

    Returns: 
    --------
    np.array - array of perturbed birth and death rates for the mutant 
    species of interest (shape: (1, 2))
    '''

    return birth_death + random.uniform(low = -0.1, high = 0.1, size = (1,2))

def pre_interaction_matrix(S: int, method = "May"): 
    '''
    Generate an interaction matrix for the relationships between species 
    using either the method described by May or Alessina-Tang

    Parameters: 
    -----------
    int S - total number of species; will generate a total of S(S-1) interaction
    parameters 
    string method - the method to generate the interaction matrix 

    Returns: 
    --------
    np.array - interaction matrix describing the relationship between the S
    species (shape: (S, S), diagonal elements are set to 0)
    '''

    # generate the interaction matrix with all entries filled 
    interaction_matrix = random.normal(loc = 0, scale = 1, size = (S, S))

    # set the diagonal elements to 0
    return np.fill_diagonal(interaction_matrix, 0)

def mutant_interaction_matrix(interaction_row: np.array, interaction_col: np.array): 
    '''
    Perturb the given interaction matrix row and column from draws 
    of a Uniform distribution on the interval [-0.1, 0.1]

    Parameters: 
    -----------
    np.array interaction_row - row of the interaction matrix for the 
    species of interest; row of the interaction matrix to be perturbed
    np.array interaction_col - column of the interaction matrix for
    the species of interest; column of the interaction matrix to be 
    perturbed

    Returns: 
    --------
    np.array (1) - array of perturbed interaction column information 
    (shape: (n-1,))
    np.array (2) - array of perturbed interaction row information 
    (shape: (n,))
    '''

    # perturb the interaction column first, dropping the last term (0)
    interaction_col = interaction_col[:-1] + \
                      random.uniform(low = -0.1, high = 0.1, 
                      size = (interaction_col.shape[0]-1))

    # perturb the interaction row
    interaction_row = interaction_row + \
                      random.uniform(low = -0.1, high = 0.1, 
                      size = (interaction_row.shape[0]))

    return interaction_row, interaction_col

def generate_interaction_matrix(S: int, method = "May"): 
    '''
    Generate an interaction matrix for the relationships between all
    species and the mutant of the species of interest using either 
    the method described by May or Alessina-Tang

    Parameters: 
    -----------
    int S - total number of species; will generate a total of S(S+1) interaction
    parameters
    string method - the method to generate the interaction matrix 

    Returns: 
    --------
    np.array - interaction matrix describing the relationship between the S
    species (shape: (S+1, S+1), diagonal elements are set to 0)
    '''

    # generate the interaction matrix without the mutant
    pre_interaction_matrix = pre_interaction_matrix(S, method = method)

    # perturb the species of interest (the last species)
    interaction_row, interaction_col = mutant_interaction_matrix(
                                        pre_interaction_matrix[-1,:],
                                        pre_interaction_matrix[:,-1])

    # stack the perturbed interactions onto the interaction matrix
    interaction_matrix = np.vstack([np.hstack([pre_interaction_matrix, 
                                               interaction_col]), 
                                    interaction_row])

    return interaction_matrix
