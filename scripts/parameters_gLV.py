
__author__ = "Ammaar A. Saeed"
__email__ = "aasaeed@college.harvard.edu"

'''
Description: 
------------
This file contains functions that generate parameters for the generalized 
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
from numpy.random import default_rng

def generate_birth_death_rates(
        S: int, 
        seed: int
    ): 
    '''
    Generate a pair of positive birth and death rates for each species
    
    Parameters: 
    -----------
    int S - total number of species; the number of birth and death rates to 
    generate
    int seed - random seed integer 

    Returns: 
    --------
    np.array - array of birth (column 1) and death (column 2) rates for each 
    species (shape: (S, 2))
    '''

    # define a rng object
    rng = default_rng(seed)

    return np.abs(rng.normal(loc = 0, scale = 1, size = (S, 2)))

def mutant_birth_death_rate(
        birth_death: np.array, 
        seed: int
    ): 
    '''
    Perturb the given birth and death rates from draws of a Uniform distribution 
    on the interval [-0.1, 0.1]

    Parameters: 
    -----------
    np.array birth_death - birth and death rates of the species of interest; 
    birth and death rates to be perturbed 
    int seed - random seed integer

    Returns: 
    --------
    np.array - array of perturbed birth and death rates for the mutant 
    species of interest (shape: (1, 2))
    '''

    # define a rng object
    rng = default_rng(seed)

    return birth_death + rng.uniform(low = -0.1, high = 0.1, size = (1,2))

def pre_interaction_matrix(
        S: int, 
        seed: int, 
        method = "May", 
        rho = None, 
        tol = None
    ): 
    '''
    Generate an interaction matrix for the relationships between species 
    using either the method described by May or Alessina-Tang

    Parameters: 
    -----------
    int S - total number of species; will generate a total of S(S-1) interaction
    parameters 
    int seed - random seed integer
    string method - the method to generate the interaction matrix: the two
    options are either (1) "May" or (2) "Alessina-Tang"
    float rho - correlation coefficient between bivariate Normal random 
    variables 

    Returns: 
    --------
    np.array - interaction matrix describing the relationship between the S
    species (shape: (S, S), diagonal elements are set to 0)
    '''

    # define a rng object
    rng = default_rng(seed)

    interaction_matrix = np.zeros((S, S))

    # using the May method
    if (method == "May"): 

        # generate the interaction matrix with all entries filled 
        interaction_matrix = rng.normal(loc = 0, scale = 1, size = (S, S))

    elif (method == "Alessina-Tang"): 

        if (rho == None): 
            raise ValueError

        # generate paired data from a bivariate normal distribution 
        pairs = rng.multivariate_normal(mean = np.array([0, 0]), 
                                        cov = np.array([1, rho, rho, 1]).reshape(2, 2),
                                        size = S * (S - 1) // 2)

        # place the paired data into the matrix 
        interaction_matrix = np.zeros((S, S))
        interaction_matrix[np.triu_indices(S, 1)] = pairs[:,0]
        interaction_matrix[np.tril_indices(S, -1)] = pairs[:,1]

    return interaction_matrix

def mutant_interaction_matrix(
        interaction_row: np.array, 
        interaction_col: np.array, 
        seed: int
    ): 
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
    int seed - random seed integer

    Returns: 
    --------
    np.array (1) - array of perturbed interaction column information 
    (shape: (n-1,))
    np.array (2) - array of perturbed interaction row information 
    (shape: (n,))
    '''

    # define a rng object
    rng = default_rng(seed)

    # perturb the interaction column first, dropping the last term (0)
    interaction_col = interaction_col + \
                      rng.uniform(low = -0.1, high = 0.1, 
                      size = (interaction_col.shape[0]))

    # perturb the interaction row
    interaction_row = interaction_row + \
                      rng.uniform(low = -0.1, high = 0.1, 
                      size = (interaction_row.shape[0]))
    interaction_row = np.hstack([interaction_row, [0]])

    return interaction_row, interaction_col

def generate_interaction_matrix(
        S: int, 
        seed: int, 
        method = "May", 
        rho = None
    ): 
    '''
    Generate an interaction matrix for the relationships between all
    species and the mutant of the species of interest using either 
    the method described by May or Alessina-Tang

    Parameters: 
    -----------
    int S - total number of species; will generate a total of S(S+1) interaction
    parameters
    int seed - random seed integer
    string method - the method to generate the interaction matrix: the two
    options are either (1) "May" or (2) "Alessina-Tang"
    float rho - correlation coefficient between bivariate Normal random 
    variables 

    Returns: 
    --------
    np.array - interaction matrix describing the relationship between the S
    species (shape: (S+1, S+1), diagonal elements are set to 0)
    '''

    interaction_matrix = None

    # generate the interaction matrix without the mutant
    interaction_matrix = pre_interaction_matrix(
                            S, 
                            method = method, 
                            seed = seed, 
                            rho = rho
                            )

    # perturb the species of interest (the last species)
    interaction_row, interaction_col = mutant_interaction_matrix(
                                        interaction_matrix[-1,:],
                                        interaction_matrix[:,-1], 
                                        seed = seed
                                        )

    # stack the perturbed interactions onto the interaction matrix
    interaction_matrix = np.vstack([np.hstack([interaction_matrix, 
                                            interaction_col[:,None]]), 
                                    interaction_row[None,:]])

    return interaction_matrix
