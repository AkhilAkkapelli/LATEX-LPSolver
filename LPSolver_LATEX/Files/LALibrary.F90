#include "Preprocessor.F90"

MODULE LALibrary
IMPLICIT NONE
!!! Linear Algebra Library


PRIVATE
PUBLIC :: DIAG, ONES, COLMULT, DOT, ENORM, GEMV, GEMM, PDSOL, TRANS, NSQRT, NLOG, ADD

CONTAINS

PURE FUNCTION DIAG(x) RESULT(D)
! Diagnoal Matrix D of vector x

UREAL, INTENT(IN) :: x(:)
UREAL             :: D(size(x), size(x))

INTEGER :: i


D=0
DO i = 1,size(x)
D(i,i) = x(i)
END DO

END FUNCTION DIAG

PURE FUNCTION ONES(n) RESULT(D)
! Diagnol Matrix D of integer dimension n

INTEGER, INTENT(IN) :: n
UREAL               :: D(n,n)

INTEGER :: i


D=0
DO i = 1,n
D(i,i) = 1
END DO

END FUNCTION ONES

PURE FUNCTION COLMULT(u,A) RESULT(uA)
! Multiply Columns of Matrix A with vector u to get Matrix uA

UREAL, INTENT(IN) :: u(:),A(:,:)
UREAL             :: uA(size(A,1),size(A,2))

INTEGER :: i


DO i = 1,size(u)
	uA(:,i) = u(i)*A(:,i)
END DO

END FUNCTION COLMULT


PURE FUNCTION ADD(v) RESULT(a)
! Add all the elements of vector v to get a scalar a

UREAL, INTENT(IN) :: v(:)
UREAL             :: a

INTEGER :: i


a = 0
DO i = 1, size(v)
  a = a + v(i)
END DO

END FUNCTION ADD


PURE FUNCTION DOT(u,v) RESULT(uTv)
! Dot Product vectors u and v t get a scalar uTv

UREAL, INTENT(IN) :: u(:),v(:)
UREAL             :: uTv


uTv = ADD(u(:)*v(:))

END FUNCTION DOT

ELEMENTAL FUNCTION NSQRT(n) RESULT(s)
! Square root s of n

UREAL, INTENT(IN) :: n
UREAL             :: s

INTEGER :: i


!IF(n<0) STOP "NSQRT ERROR: Invalid Input"

s = n/2

DO i = 1,realKind
  s = (s + n/s)/2
END DO

END FUNCTION NSQRT

ELEMENTAL FUNCTION NLOG(n) RESULT(l)
! Natural Logarithm l of n

UREAL, INTENT(IN) :: n
UREAL             :: l


l = 2*ATANH((n-1)/(n+1))

END FUNCTION NLOG

PURE FUNCTION ENORM(v) RESULT(n)
! Euclidean Norm scalar n of vector v

UREAL, INTENT(IN) :: v(:)
UREAL             :: n

n = NSQRT(DOT(v,v))

END FUNCTION ENORM

PURE FUNCTION CholeskyDecomposition(A) RESULT(L)
! Lower Triangular Matrix L by Cholesky Decomposition of Matrix A

UREAL, INTENT(IN) :: A(:,:)
UREAL             :: L(size(A,1),size(A,2))

INTEGER :: i
UREAL   :: B(size(A,1), size(A,2)), summ


!IF(size(A,1) /= size(A,2)) STOP "CHOLESKY ERROR: Invalid size"
B = A
L = 0

DO i = 1,size(B,1) 
summ = B(i,i) - DOT(B(i,:i-1),B(i,:i-1))
!IF(summ <= 0.) STOP "CHOLESKY ERROR: Invalid Matrix Input"
L(i,i) = NSQRT(summ)
B(i+1:,i)=(B(i,i+1:)-GEMV(B(i+1:,:i-1),B(i,:i-1)))/L(i,i)
L(i+1:,i) = B(i+1:,i)
END DO

END FUNCTION CholeskyDecomposition

PURE FUNCTION ForSubstitution(L, u) RESULT(v) 
! vector v by Forward Substitution of Lower Triangular Matrix L and vector u

UREAL, INTENT(IN) :: L(:,:), u(:)
UREAL             :: v(size(u,1))

INTEGER :: i


!IF(size(L,2) /= size(u)) STOP "BackSubstitution ERROR: Invalid size input"

DO i  = 1,size(u,1)
v(i) = (u(i) - ADD(L(i,:i-1)*v(:i-1)))/L(i,i)
END DO

END FUNCTION ForSubstitution

PURE FUNCTION BackSubstitution(U, a) RESULT(v) 
! vector v by Forward Backward of Upper Triangular Matrix U and vector a

UREAL, INTENT(IN) :: U(:,:), a(:)
UREAL             :: v(size(a,1))

INTEGER :: i


!IF(size(U,2) /= size(b)) STOP "ForSubstitution ERROR: Invalid size input"

DO i  = size(a,1),1,-1
v(i) = (a(i) - SUM(U(i,i+1:)*v(i+1:)))/U(i,i)
END DO

END FUNCTION BackSubstitution


PURE FUNCTION GEMV(M,v) RESULT(Mv)
! vector Mv by Matrix Vector Multiplication of Matrix M and vector v

UREAL, INTENT(IN) :: M(:,:),v(:)
UREAL             :: Mv(size(M,1))

INTEGER :: i


!IF(size(M,2) /= size(v)) STOP "MatVecMult ERROR: Invalid size input"

DO i = 1,size(M,1)
	Mv(i) = DOT(M(i,:),v(:))
END DO

END FUNCTION GEMV

PURE FUNCTION GEMM(A,B) RESULT(AB)
! Matrix AB by Matrix Vector Multiplication of Matrix A and Matrix B

UREAL, INTENT(IN) :: A(:,:),B(:,:)
UREAL             :: AB(size(A,1),size(B,2))

INTEGER :: i, j


DO i = 1,size(A,1)
	DO j = 1,size(B,2)
		AB(i,j) = DOT(A(i,:),B(:,j))
	END DO
END DO

END FUNCTION GEMM

PURE FUNCTION PDSOL(A, y) RESULT(x)
!vector x by solving System Ax = y of Positive Definite Symmetric Matrix A and vector y

UREAL, INTENT(IN) :: A(:,:), y(:)
UREAL             :: x(size(A,2))

UREAL :: L(size(A,1), size(A,1)), u(size(A,1))


!IF(size(A, 1) /= size(A, 2) .OR. size(A,1) /= size(y)) STOP "PDLSS ERROR: Invalid Input"

L = CholeskyDecomposition(A)

u = ForSubstitution(L, y)

x = BackSubstitution(TRANS(L), u)

END FUNCTION PDSOL

PURE FUNCTION TRANS(A) RESULT(AT)
! Transpose AT of Matrix A

UREAL, INTENT(IN) :: A(:,:)
UREAL             :: AT(size(A,2),size(A,1))

INTEGER :: i, j


FORALL(i = 1:size(A,1) , j = 1:size(A,2))
	AT(j,i) = A(i,j)
END FORALL

END FUNCTION TRANS

END MODULE LALibrary
