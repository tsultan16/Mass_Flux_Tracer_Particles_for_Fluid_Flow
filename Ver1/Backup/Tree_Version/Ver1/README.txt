# Mass-Flux/"Monte Carlo" Tracers

# Algorithm for tracer advection based on Genel et al., 2013, MNRAS, 435, 1426  

#Core Modules:   1) Solver Module
                 2) Data Structure Module


# Some Notes about the Tree Data Structure: Inside each cell, tracer particles are divided among clusters.
                                            Each cluster is a branch/node of the cell's "tree". Balance can be
									        maintained (or imbalances can be prevented from becoming worse) by
											always placing new tracers in the branch with fewest tracers.


									        When picking outgoing tracers (uniform) randomly from an imbalanced tree,
											there will be bias since tracers resiing in branches with higher population
											have higher probabilty of getting picked. To avoid this bias, need to pick tracers
											using a non-uniform random distribution (more weight on branches with lower population)


# Some Notes about Performance:    Most of the computational cost comes from 
                                   the "cell_transfer" subroutine, which is called
                                   O(N) times for a cell containing N tracers.
                                   The cpu_time for each cell_face cycle comes out to
                                   about O(N^2).
        
