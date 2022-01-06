#include "Preprocessor.F90"

MODULE LPTools
USE AUGOperator, ONLY : OPERATOR(.VAUG.)
USE ProbabilityGeneration, ONLY : ProbGen
USE LALibrary, ONLY : ONES, DOT, ColMult, GEMV, ENORM, TRANS, NLOG, ADD, NSQRT
IMPLICIT NONE

PRIVATE
PUBLIC :: Dual, StdToCan, ProjectiveTransform, InvProjTransform, zerostop, potentialstop, optimumstop, zeroratio, minratio, potentialratio

CONTAINS


UREAL PURE FUNCTION Potential(c,x) RESULT(f)

UREAL, INTENT(IN) :: c(column+1) , x(column+1)


f = ADD(NLOG(DOT(c,x)/x))

END FUNCTION Potential

LOGICAL PURE FUNCTION zerostop(n, x,xp, iter, c) RESULT(stp)

INTEGER, INTENT(IN) :: iter, n
UREAL, INTENT(IN) :: x(n), xp(n), c(n)

INTEGER, PARAMETER :: iterlimit = 10000


stp = .FALSE.

IF(iter >= iterlimit .OR. x(n-1) < real(1.Q-100,realKind) .OR. x(n-1) > xp(n-1)) stp = .TRUE.

END FUNCTION zerostop

LOGICAL PURE FUNCTION potentialstop(n, x,xp, iter, c) RESULT(stp)

INTEGER, INTENT(IN) :: iter, n
UREAL, INTENT(IN) :: x(n), xp(n), c(n)

UREAL, PARAMETER :: delta= -one/8
INTEGER, PARAMETER :: iterlimit = 10000
UREAL ::  f, fp


stp = .FALSE.
fp = Potential(c,xp)
f =  Potential(c,x)

IF(iter >= iterlimit .OR. f - fp > delta) stp = .TRUE.

END FUNCTION potentialstop

LOGICAL PURE FUNCTION optimumstop(n, x,xp, iter, c) RESULT(stp)

INTEGER, INTENT(IN) :: iter, n
UREAL, INTENT(IN) :: x(n), xp(n), c(n)

INTEGER, PARAMETER :: iterlimit = 10000
UREAL :: obj, objp

stp = .FALSE. 
obj = DOT(c,x)
objp = DOT(c,xp)

IF(iter >= iterlimit .OR. objp - obj < real(1.Q-100,realKind) ) stp = .TRUE.

END FUNCTION optimumstop

UREAL PURE FUNCTION potentialratio(n, cunit) RESULT(alpha)

INTEGER, INTENT(IN) :: n
UREAL, INTENT(IN) :: cunit(n)


alpha = 1/(4*NSQRT(ukind(n)*(ukind(n)-1)))

END FUNCTION potentialratio

UREAL PURE FUNCTION zeroratio(n, cunit) RESULT(alpha)

INTEGER, INTENT(IN) :: n
UREAL, INTENT(IN) :: cunit(n)

UREAL, PARAMETER :: beta = ukind(1.Q-1)
UREAL :: a
INTEGER :: idx 


alpha = one/(n*cunit(n-1))
DO idx=1,n
  IF(cunit(idx)<0 .OR. idx == n-1) CYCLE
  a = (ukind(1.Q0)-beta)/(n*cunit(idx))
  IF(alpha>a) alpha = a
END DO

END FUNCTION zeroratio

UREAL PURE FUNCTION minratio(n, cunit) RESULT(alpha)

INTEGER, INTENT(IN) :: n
UREAL, INTENT(IN) :: cunit(n)

UREAL, PARAMETER :: beta = ukind(1.Q-1)


alpha = (one-beta)/(n*maxval(cunit))

END FUNCTION minratio


PURE SUBROUTINE Dual(A, b, c, Adual, bdual, cdual)

UREAL, INTENT(IN) :: A(:,:), b(:), c(:)

UREAL, INTENT(OUT) :: Adual(size(A,2),size(A,1)), bdual(size(c)), cdual(size(b))


cdual = -b

Adual = -TRANS(A)

bdual = -c

END SUBROUTINE Dual

PURE SUBROUTINE StdtoCan(Astd,bstd,cstd ,Acan,ccan, a0)

UREAL, INTENT(IN) :: Astd(:,:), bstd(:), cstd(:)
UREAL, INTENT(OUT) :: Acan(size(Astd,1)+size(Astd,2)+1,2*(size(Astd,1)+size(Astd,2)+1)), &
				 ccan(2*(size(Astd,1)+size(Astd,2)+1)), a0(2*(size(cstd)+size(bstd))+1)

UREAL :: A(size(Astd,1)+size(Astd,2)+1,2*(size(Astd,1)+size(Astd,2))+1), &
			b(size(Astd,1)+size(Astd,2)+1),c(2*(size(Astd,1)+size(Astd,2))+1)
UREAL :: x0(size(Astd,2)), y0(size(Astd,1)), u0(size(Astd,1)), v0(size(Astd,2)), lambda0

INTEGER :: m, n, i,j


m = size(Astd,1); n = size(Astd,2)

x0 = 1
y0 = 1
u0 = 1 
v0 = 1
lambda0 = 1

a0 = x0 .VAUG. y0 .VAUG. u0 .VAUG. v0 .VAUG. lambda0

A = 0

A(1:m,1:n) = Astd
A(1:m,n+1:n+m) = -Ones(m)

A(m+1:m+n,m+n+1:2*m+n) = TRANS(Astd)
A(m+1:m+n,2*m+n+1:2*(m+n)) = Ones(n)

A(m+n+1,1:n) = cstd
A(m+n+1,m+n+1:2*m+n) = -bstd

A(1:m,2*(m+n)+1) =  bstd - GEMV(Astd,x0) + y0 
A(m+1:m+n,2*(m+n)+1) = cstd - GEMV(TRANS(Astd),u0) - v0 
A(m+n+1,2*(m+n)+1) = -DOT(cstd,x0) + DOT(bstd,u0) 


b(:m) = bstd 
b(m+1:m+n) = cstd
b(m+n+1) = 0


c = 0
c(2*m+2*n+1) = 1

CALL ProjectiveTransform(A,b,c, Acan,ccan, a0)

END SUBROUTINE StdToCan

PURE SUBROUTINE ProjectiveTransform(A,b,c, Acan,ccan, a0)

UREAL, INTENT(IN) :: A(:,:),b(:),c(:), a0(:)

UREAL, INTENT(OUT) :: Acan(size(A,1),size(A,2)+1), ccan(size(c)+1)


Acan(:,:size(A,2)) = ColMult(a0,A)
Acan(:,size(A,2)+1) = -b

ccan(:size(c)) = a0*c
ccan(size(c)+1) = 0

END SUBROUTINE ProjectiveTransform

PURE FUNCTION InvProjTransform(xcan,x0) RESULT(x)
	
UREAL, INTENT(IN) :: xcan(:), x0(:)

UREAL             :: x(size(xcan)-1)


x = (x0*xcan(:size(x0)))/xcan(size(xcan))

END FUNCTION InvProjTransform

END MODULE LPTools
