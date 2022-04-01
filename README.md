# lotka-volterra-evolution
This repository contains code used to simulate the generalized Lotka-Volterra equations for a system that undergoes evolution.

# Goals 
Start with a n=2 system for two initial species with some interaction matrix (neutral, competitive, or cooperation) and then introduce a mutation into one of the two species
Slightly perturb the interaction row of the mutation and then simulate from here/solve analytically using the Lotka-Volterra model 
Simulate this using a stochastic framework such that there are N initial species of each type and then use birth-death or some other model for how evolution will occur

# Questions
How do we choose an interaction matrix (consider a Normal [2 params] or Uniform [2 params; -0.1 and 0.1])
How do we perturb the interaction row once the mutation has occurred. 

# Computational 
Generate interaction matrix 
Perturbation 
