PROGRAM xorshift32
IMPLICIT NONE

INTEGER (KIND = 4):: s ! 4-byte signed integers
INTEGER (KIND = 4), PARAMETER :: imin = -2147483648
INTEGER (KIND = 4), PARAMETER :: imax = 2147483647
REAL, PARAMETER :: ran_inv_32 = 2.32830643653869628906E-10 ! ran_inv_32 = 2^-32 in floating point
INTEGER, PARAMETER :: N = 100000000
REAL :: p, t1, t2, t

INTEGER :: i

! set a non-zero seed value
s = 164972


t = 0.0

DO i = 1,N
   CALL CPU_TIME(t1)
   p = rand_num()
   CALL CPU_TIME(t2)
   t = t + t2 - t1
END DO
 

PRINT*,'XOR SHIFT (32 bit):'
PRINT*,'N, Time(sec) = ',N,t


CONTAINS


! returns random number between 0 and 1
FUNCTION rand_num() RESULT(r)

REAL :: r
INTEGER (KIND = 4):: x


! 32-bit xor shift 
x = IEOR(s,ISHFT(s,13))
x = IEOR(x,ISHFT(x,-17))
s = IEOR(x,ISHFT(x,5))

r=0.5+s*ran_inv_32

END FUNCTION rand_num


END PROGRAM xorshift32
