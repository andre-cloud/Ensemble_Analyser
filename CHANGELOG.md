# Version 0.3.0

## Major updates:

1. Normalize the electronic graph only in the window that has to be considered
2. Created a stand-alone script to generate the graphs even after all the calculation ended

## Further updates:

1. Clustering based on dihedral rather than on absolute position
2. Clustering method
3. IR graphs
4. Protocol keyword "No cut", to not prune the ensemble

# Version 0.2.0

## Major updates:

1. New optimizer logic, based on BFGS algorithm implemented in ASE. This action is for scalability purpose
    - New parameters for protocol definition: 
        - FMAX = max force value to consider the optimization converged
        - CONSTRAINS = list of indexes to achieve a constrained optimization
        - MAXSTEP = max displacement during the optimization

2. Each conformer has been assigned to a color, in order to be recognized in the PCA graph
3. PCA now shows a color-map interpolation of the energy surface 

## Further modification: 

1. Normalize the electronic graph only in the window that has to be considered
2. Write a complete guide
3. Benchmark the best values for the MAXSTEP and FMAX parameter.
