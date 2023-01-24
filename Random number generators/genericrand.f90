PROGRAM randGeneric
IMPLICIT NONE

INTEGER, PARAMETER :: N = 100000000

REAL :: p,t1,t2,t

INTEGER :: i

t = 0.0

DO i = 1,N
   CALL CPU_TIME(t1)
   CALL RANDOM_NUMBER(p)
   CALL CPU_TIME(t2)
   t = t + t2 - t1
END DO
 
PRINT*,'Fortran Intrinsic Random:'
PRINT*,'N, Time(sec) = ',N,t



END PROGRAM randGeneric
