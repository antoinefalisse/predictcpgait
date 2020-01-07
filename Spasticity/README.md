# Spasticity

This repository contains code to derive personalized spasticity models. This code is inspired from [Falisse et al. (2018)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0208811).

Deriving personalized spasticity models consists of different steps.

1. NormalizeIPSAEMG
    - This folder contains code to normalize the EMG data from the IPSA measurements.
2. ForwardSimulation
    - This folder contains code to generate EMG-driven forward simulations.
3. SpasticGainEstimation
    - This folder contains code to estimate personalized feedback gains based on passive movements (IPSA measurements).
4. SpasticityGait
    - This folder contains code to predict the spastic response during gait.
    
Each step depends on results from the previous one, please keep in mind the order. Furthermore, the spasticity models depend on the muscle-tendon parameters(see [ParameterEstimation](https://github.com/antoinefalisse/predictcpgait/tree/master/ParameterEstimation)).