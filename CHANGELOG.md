# Version 0.2.0

## Major updates:

1. New optimizer logic, based on BFGS algorithm implemented in ASE. This action is for scalability purpose
    - New parameters for protocol definition: 
        - FMAX = max force value to consider the optimization converged
        - CONSTRAINS = list of indexes to achieve a constrained optimization
        - MAXSTEP = max displacement during the optimization

2. Each conformer has been assigned to a color, in order to be recognized in the PCA graph

## Further modification: 

1. Normalize the electronic graph only in the window that has to be considered
2. Consider the possibility to turn on a color-map visualization in the PCA
3. Write a complete guide
4. Benchmark the best values for the MAXSTEP and FMAX parameter.
