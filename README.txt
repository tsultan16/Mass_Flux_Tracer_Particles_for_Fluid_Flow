# Mass-Flux/"Monte Carlo" Tracers

# Algorithm for tracer advection based on Genel et al., 2013, MNRAS, 435, 1426  
# Algorithm for balanced tree based on Knuth, Art of C. P. Vol. 3

#Core Modules:   1) Solver Module
                 2) Data Structure Module


# Some Notes about the Tree Data Structure: Inside each cell, tracers are assembled in a binary search tree
                                            data structure. Height balance is maintained using the AVL tree algorithm.
                                            Tracers are sorted according to their "rank" attribute.
                                             

# Some Notes about Performance: The cost of inserting or removing a tracer from it's tree is O[log_2(N)]
                                (much faster than the straight up linked list version).
        


# Implemented 1D upto 3D

# Need to make a new subroutine for "creating" and "destroying" tracer particles (i.e. allocating and deallocating memory)
