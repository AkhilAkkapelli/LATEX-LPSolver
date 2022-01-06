#include "Preprocessor.F90"

MODULE ProjScalingAlgo
USE AUGOperator, ONLY : OPERATOR(.VAUG.)
USE LALibrary, ONLY : DOT, ENORM, GEMV, GEMM, TRANS, PDSOL, COLMULT
IMPLICIT NONE


PUBLIC :: ProjectiveScale
PRIVATE

CONTAINS


PURE FUNCTION ProjectiveScale(A,c, stopcond, ratiotest) RESULT(xcan)
! Obtain vector xcan by minimizing along Objective vector c of system Axcan = 0 using stoping condition stopcond and ratio test ratiotest functions

UREAL, INTENT(IN) :: A(:,:), c(:)
UREAL             :: xcan(size(A,2))

UREAL :: xp(size(A,2)), x(size(A,2))
INTEGER :: n, iter

INTERFACE 
  
  UREAL PURE FUNCTION ratiotest(n, cunit) RESULT(alpha)
    INTEGER, INTENT(IN) :: n
    UREAL,   INTENT(IN) :: cunit(n)
	END FUNCTION ratiotest
	
	LOGICAL PURE FUNCTION stopcond(n, x,xp, iter, c) RESULT(stp)
  INTEGER, INTENT(IN) :: iter, n
  UREAL, INTENT(IN) :: x(n), xp(n), c(n)
  END FUNCTION stopcond
  
END INTERFACE


!IF(size(A,2) /= size(c)) STOP "Algorithm ERROR: Wrong size for input argument"

n = size(A,2)
xp = one/n

!IF( ANY(GEMV(A,xp) >= ukind(1.Q-10)) ) STOP "Center of Simplex doesn't lie in Null space of the Input Matrix:"

iter = 1
DO 

x = Optimize(xp)
IF(stopcond(n, x,xp, iter, c)) EXIT

!IF( ANY(ISNAN(x))) THEN
!PRINT*,  "WARNING: NAN occurred in the Solution"
!x= xp
!EXIT 
!END IF

xp = x
iter = iter + 1

END DO	

xcan = x

CONTAINS 

PURE FUNCTION Optimize(xp) RESULT(x)

UREAL, INTENT(IN) :: xp(:)
UREAL             :: x(size(xp))

UREAL :: Ad(size(A,1),size(xp)), B(size(A,1)+1,size(xp)), &
	      	v(size(A,1)+1), cp(size(c)), cunit(size(c)), alpha
UREAL ::  e(size(xp)) = one 
UREAL :: x0(size(xp)) = e/n


Ad = COLMULT(xp,A)
B = Ad .VAUG. e
v = PDSOL(GEMM(B,TRANS(B)), GEMV(B,xp*c))

cp = xp*c - GEMV(TRANS(B),v)
cunit = cp/ENORM(cp)

alpha = ratiotest(n, cunit)

x = x0 - alpha*cunit
x = x*xp/DOT(x,xp)

END FUNCTION Optimize

END FUNCTION ProjectiveScale


END MODULE ProjScalingAlgo
