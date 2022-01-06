#include "Preprocessor.F90"

MODULE LinearProblemSolver
USE ProbabilityGeneration, ONLY : ProbGen 
USE ProjScalingAlgo, ONLY : ProjectiveScale
USE LPTools, ONLY : StdToCan, ProjectiveTransform, InvProjTransform, Dual, zerostop, potentialstop, optimumstop, zeroratio, minratio, potentialratio
USE LALibrary, ONLY : DOT, TRANS, GEMV, ENORM, PDSOL, GEMM
USE AUGOperator, ONLY : OPERATOR(.VAUG.), OPERATOR(.HAUG.)
IMPLICIT NONE

PRIVATE
PUBLIC :: LPType

CONTAINS

SUBROUTINE Canonical

UREAL :: Acan(row,column), ccan(column), xopt(column)
INTEGER :: i, j


OPEN(unit = 10, file= './1.IO/Acan.txt')
  DO i=1,row
    READ(10, *) (Acan(i,j) ,j=1,column)
  END DO
CLOSE(10) 

OPEN(unit = 30, file= './1.IO/ccan.txt')
  DO j=1,column
    READ(30, *) ccan(j)
  END DO
CLOSE(30)

!CALL ProbGen(c)

xopt = ProjectiveScale(Acan, ccan, potentialstop, potentialratio)

END SUBROUTINE Canonical

SUBROUTINE Equality

UREAL :: Aeq(row, column), beq(row), ceq(column), x0(column), xeq(column), eps, &
			Acan(size(Aeq,1),size(Aeq,2)+1), ccan(size(ceq)+1), xcan(size(Aeq,2)+1), &
			Acan1(size(Aeq,1),size(Aeq,2)+2), ccan1(size(ceq)+2), xcan1(size(Aeq,2)+2), &
				A(row,column+1), b(row), c(column+1), x01(column+1), xopt(column), xopt1(column+1)
INTEGER :: i, j


OPEN(unit = 10, file= './1.IO/Aeq.txt')
DO i=1,row
READ(10, *) (Aeq(i,j) ,j=1,column)
END DO
CLOSE(10)  

OPEN(unit = 20, file= './1.IO/beq.txt')
DO i=1,row
     READ(20, *) beq(i)
END DO
CLOSE(20)

OPEN(unit = 30, file= './1.IO/ceq.txt')
DO j=1,column
     READ(30, *) ceq(j)
END DO
CLOSE(30) 

x0 = one

x01 = x0 .VAUG. one

b = beq + 0*(beq - GEMV(Aeq,x0))

A(:,:column) = Aeq
A(:,column+1) = 1*(beq - GEMV(Aeq,x0))

c(:column) = 0
c(column+1) = one

CALL ProjectiveTransform(A,b,c, Acan1,ccan1, x01)

!ccan = ccan - MAX(MAXVAL(ccan),0.)

xcan1 = ProjectiveScale(Acan1, ccan1, optimumstop, potentialratio)

xopt1 = InvProjTransform(xcan1, x01)

x0 = xopt1(:column) 
print*, 'x0', xopt1
IF(x0(column+1) < 1.Q0-10) THEN

CALL ProjectiveTransform(Aeq,beq,ceq, Acan,ccan, x0)

xcan = ProjectiveScale(Acan, ccan, potentialstop, potentialratio)

xopt = InvProjTransform(xcan, x0)
print*, 'xopt', xopt
END IF

END SUBROUTINE Equality

SUBROUTINE Standard

UREAL :: Astd(row,column), bstd(row), cstd(column), xopt(column), &
		      Acan(size(Astd,1)+size(Astd,2)+1,2*(size(Astd,1)+size(Astd,2)+1)), ccan(2*(size(Astd,1)+size(Astd,2)+1)), &
		       xcan(2*(size(Astd,1)+size(Astd,2)+1)), xopt1(2*(size(Astd,1)+size(Astd,2))+1), a0(2*(size(Astd,1)+size(Astd,2))+1)
INTEGER :: i, j


OPEN(unit = 10, file= './1.IO/Astd.txt')
DO i=1,row
READ(10, *) ( Astd(i,j) ,j=1,column)
END DO
CLOSE(10)

OPEN(unit = 40, file= './1.IO/bstd.txt')
DO i=1,row
     READ(40, *) bstd(i)
END DO
CLOSE(40)

OPEN(unit = 30, file= './1.IO/cstd.txt')
DO j=1,column
     READ(30, *) cstd(j)
END DO
CLOSE(30)

CALL StdtoCan(Astd,bstd,cstd, Acan,ccan, a0)

xcan = ProjectiveScale(Acan,ccan, potentialstop, potentialratio)

xopt1 = InvProjTransform(xcan, a0)

xopt = xopt1(:column)

END SUBROUTINE Standard

SUBROUTINE LeastNegative

UREAL :: Aln(row, column), bln(row), cln(column), y0(row), x0(column), xln(column), eps, e(column), &
		      Aaug(row,column+1), baug(row), caug(column+1), x0aug(column+1), xaug(column+1), &
		    	Aaugcan(size(Aln,1),size(Aln,2)+2), caugcan(size(cln)+2), xaugcan(size(Aln,2)+2), &
		  		Acan(size(Aln,1),size(Aln,2)+1), ccan(size(cln)+1), xcan(size(Aln,2)+1)	
			
INTEGER :: i, j


OPEN(unit = 10, file= './1.IO/Aln.txt')
DO i=1,row
READ(10, *) (Aln(i,j) ,j=1,column)
END DO
CLOSE(10)  

OPEN(unit = 20, file= './1.IO/bln.txt')
DO i=1,row
     READ(20, *) bln(i)
END DO
CLOSE(20)

OPEN(unit = 30, file= './1.IO/cln.txt')
DO j=1,column+1
     READ(30, *) cln(j)
END DO
CLOSE(30)

e=one

y0 = PDSOL(GEMM(Aln,TRANS(Aln)), bln)

x0 = GEMV(TRANS(Aln), y0)

IF(ANY(x0 < 0)) THEN

eps = -2*minval(x0)

Aaug = Aln .HAUG. -GEMV(Aln,e)

baug = bln 

x0aug = (x0 + eps) .VAUG. eps

caug = 0
caug(column+1) = 1

CALL ProjectiveTransform(Aaug,baug,caug, Aaugcan, caugcan, x0aug)

xaugcan = ProjectiveScale(Aaugcan, caugcan, potentialstop, potentialratio)

xaug = InvProjTransform(xaugcan, x0aug)

x0 = xaug(:column)

END IF

IF(xaug(column+1) < 1.0Q-30) THEN

CALL ProjectiveTransform(Aln,bln,cln, Acan,ccan, x0)

ccan = ccan - MAX(MAXVAL(ccan),0.)

xcan = ProjectiveScale(Acan, ccan, potentialstop, potentialratio)

xln = InvProjTransform(xcan, x0)

ELSE

eps = xaug(column +1)

END IF

END SUBROUTINE LeastNegative

END MODULE LinearProblemSolver
