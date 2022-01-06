#include "Preprocessor.F90"

#ifdef ProbDist
MODULE LPPGeneration
USE ProbabilityDistribution, ONLY : ProbDist
IMPLICIT NONE
! Generate LP Problem from Probability Distribution

PRIVATE
PUBLIC :: LPGen

CONTAINS 

SUBROUTINE LPGen
! Generate Linear Programming Problem

UREAL :: A(row,column), b(row), c(column), x0(column)
INTEGER :: seed, count
INTEGER :: i, j


CALL SYSTEM_CLOCK(count)
seed = MOD(count,10000)

OPEN(unit = 10, file= './1.IO/A.txt')
DO i = 1, row
DO j = 1, column
	A(i,j) = ProbDist(seed)
	WRITE(10, *) A(i,j)
END DO
WRITE(10, *) ""
END DO
CLOSE(10) 

OPEN(unit = 20, file= './1.IO/x0.txt')
DO i = 1, column
	x0(i) = ProbDist(seed)
	WRITE(20, *) x0(i)
END DO
CLOSE(20)

OPEN(unit = 30, file= './1.IO/c.txt')
DO i = 1, column
	c(i) = ProbDist(seed)
	WRITE(30, *) c(i)
END DO
CLOSE(30)

OPEN(unit = 40, file= './1.IO/b.txt')
DO i = 1, row
	b(i) = ProbDist(seed)
	WRITE(40, *) b(i)
END DO
CLOSE(40)

END SUBROUTINE LPGen

END MODULE LPPGeneration
#endif
