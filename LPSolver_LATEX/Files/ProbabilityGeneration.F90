#include "Preprocessor.F90"

#ifdef ProbDist
MODULE ProbabilityGeneration
USE ProbabilityDistribution, ONLY : ProbDist
IMPLICIT NONE
! Generate Probability Distribution vector/Matrix

PRIVATE
PUBLIC :: ProbGen

INTERFACE ProbGen

	MODULE PROCEDURE ProbVecGen
	MODULE PROCEDURE ProbMatGen

END INTERFACE

CONTAINS 

SUBROUTINE ProbVecGen(v, trail)
!Generate Probablity Distribution Vector v integer trail

UREAL, INTENT(OUT) :: v(:)
INTEGER, OPTIONAL, INTENT(IN)  :: trail

INTEGER :: i, j, seed, count


IF(PRESENT(trail)) THEN
  seed = 12345
  DO j = 1, trail
    DO i = 1, size(v)
    	v(i) = ProbDist(seed)
    END DO	
  END DO
ELSE
  CALL SYSTEM_CLOCK(count)
  seed = MOD(count,10000)
END IF  

DO i = 1, size(v)
	v(i) = ProbDist(seed)
END DO

END SUBROUTINE ProbVecGen

SUBROUTINE ProbMatGen(M, trail)
!Generate Probablity Distribution Matrix M with integer trail

UREAL, INTENT(OUT)             :: M(:,:)
INTEGER, OPTIONAL, INTENT(IN)  :: trail

INTEGER :: i, j, k, seed, count


IF(PRESENT(trail)) THEN
  seed = 12345
  DO k = 1, trail
    DO i = 1, size(M,1)
      DO j = 1, size(M,2)
      	M(i,j) = ProbDist(seed)
      END DO
    END DO
  END DO
ELSE
  CALL SYSTEM_CLOCK(count)
  seed = MOD(count,10000)
END IF  

DO i = 1, size(M,1)
  DO j = 1, size(M,2)
  	M(i,j) = ProbDist(seed)
  END DO
END DO

END SUBROUTINE ProbMatGen

END MODULE ProbabilityGeneration
#endif
