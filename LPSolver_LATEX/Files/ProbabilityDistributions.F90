#include "Preprocessor.F90"

#ifdef ProbDist
MODULE ProbabilityDistribution
IMPLICIT NONE

PRIVATE
PUBLIC :: ProbDist

CONTAINS

FUNCTION Uniform(seed) RESULT(u)
! Uniform Distribution scalar u from integer seed

INTEGER :: seed
INTEGER, PARAMETER :: IA=16807,IM=2147483647,IQ=127773,IR=2836
UREAL, SAVE :: am
INTEGER, SAVE :: ix=-1,iy=-1,k

UREAL :: u


! Initialise
IF (seed <= 0 .OR. iy < 0) THEN
am=nearest(1.0,-1.0)/IM
iy=ior(ieor(888889999,abs(seed)),1)
ix=ieor(777755555,abs(seed))
seed=abs(seed)+1
END IF

!!! Marsaglia shift sequence with period 2^32 -1
ix=ieor(ix,ishft(ix,13))
ix=ieor(ix,ishft(ix,-17))
ix=ieor(ix,ishft(ix,5))

! Park-Miller sequence by Schrage’s method with period2^31−2
k=iy/IQ
iy=IA*(iy-k*IQ)-IR*k
if (iy < 0) iy=iy+IM

! Combine two Generators
u=am*ior(iand(IM,ieor(ix,iy)),1)

#ifdef lowerbound 
#ifdef upperbound
u = (upperbound-lowerbound)*u + lowerbound
#endif
#endif

END FUNCTION Uniform

FUNCTION Normal(seed) RESULT(n)
! Normal Distribution scalar n from integer seed

INTEGER :: seed
UREAL, PARAMETER :: PI = 4.Q0*DATAN(1.D0)
UREAL :: u1,u2

UREAL :: n


u1 =  Uniform(seed)
u2 =  Uniform(seed)

n = sqrt(-2*log(u1))*cos(2*PI*u2)

#ifdef mean
#ifdef variance
n = variance*n + mean
#endif
#endif

END FUNCTION Normal

END MODULE ProbabilityDistribution
#endif
