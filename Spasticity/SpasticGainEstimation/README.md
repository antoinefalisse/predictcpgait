# SpasticGainEstimation

This repository contains code to estimate personalized feedback gains based on passive movements (IPSA measurements). It consists of two main script:

1. SpasticGainEstimation_main.m
    - This script estimates feedback gains to reproduce data from passive movements (IPSA measurements).
2. SpasticGainEstimation_observeResults.m
    - This script allows observing the results of the feedback gain estimation.
    
The two other scripts are optional:

1. adjustEMGOnset.m
    - This script allows adjusting the EMG onset identified automatically in the IPSA data.
2. selectOptimizationRange.m
    - This script allows adjusting the optimization time window used when estimating the feedback gains.