! Fortran version of Taus88 C-code from L’ECUYER(1996)

PROGRAM taus88

IMPLICIT NONE

INTEGER, PARAMETER :: N = 100000000

! 32 bit integers
INTEGER(KIND=4) :: s1, s2, s3, b
REAL :: dat, t1,t2,t

INTEGER :: i


! set non-zero seed values
s1 = 1234
s2 = 5678
s3 = 9123



t = 0.0

DO i= 1, N
        CALL CPU_TIME(t1)
        dat=rand_num()
        CALL CPU_TIME(t2)
        t = t + t2 - t1    
END DO


PRINT*,'TAUS88:'
PRINT*,'N,Time(sec)=',N,T

CLOSE(UNIT = 10)


CONTAINS


! returns random number between 0 and 1
FUNCTION rand_num() RESULT(r)

    REAL :: r


    b = ISHFT(IEOR(ISHFT(s1,13),s1),-19)    
    s1 = IEOR(ISHFT(IAND(s1,4294967294),12),b)
    b = ISHFT(IEOR(ISHFT(s2,2),s2),-25)
    s2 = IEOR(ISHFT(IAND(s2,4294967288),4),b)
    b = ISHFT(IEOR(ISHFT(s3,3),s3),-11)
    s3 = IEOR(ISHFT(IAND(s3,4294967280),17),b)

    r = 0.5+IEOR(s1,IEOR(s2,s3))*2.3283064365E-10
 


END FUNCTION rand_num

END PROGRAM taus88
