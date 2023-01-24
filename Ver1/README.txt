# Mass-Flux/"Monte Carlo" Tracers

# Algorithm for tracer advection based on Genel et al., 2013, MNRAS, 435, 1426  

#Core Modules:   1) Solver Module
                 2) Data Structure Module


# Some Notes about the Tree Data Structure: Inside each cell, tracer particles are divided among clusters.
                                            Each cluster is a branch/node of the cell's "tree". (Weight) Balance can be
									        maintained (or imbalances can be prevented from becoming worse) by
											always placing new tracers in the branch with fewest tracers.


									        When picking outgoing tracers (uniform) randomly from an imbalanced tree,
											there will be bias since tracers resiing in branches with higher population
											have higher probabilty of getting picked. To avoid this bias, need to pick tracers
											using a non-uniform random distribution (more weight on branches with lower population)


# Some Notes about Performance: Suppose each tree in a cell has p branches. The tree is weight balanced,
                                so there are N/p tracers in each branch. The number of steps required to 
                                find a given tracer is given by:

                                s = a*(N/p) + b*p  , where  a,b are constants and 0< a <= p, 0 < b <= N/p       

                                Then to minimize s with respect to p:

                                ds/dp = -a*(N/p^2) + b = 0 => p_min = (a*N/b)^1/2 

                                => p_min = Order(N^1/2) => s = Order(N^1/2)
                                
        
