# This is an augmented version of Ver2. It is parallelzied using MPI domain decomposition.
  It includes new pack and unpack subroutines for transferring tracers between MPI ranks.

# Need to try:

1) Make the boundary cell data structure a stack/linked list rather than a tree since all tracers 
   contained inside are going to be copied over to the neighbor rank. This should make things simpler. 