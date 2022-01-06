FUNCTION Karmarkar(A,c) RESULT(x_can)

REAL, INTENT(IN) :: A(:,:), c(:)
REAL :: x_can(size(c))
REAL :: x_p(size(A,2)), x(size(A,2)), f_p, f, delta, eps
INTEGER :: n


IF(size(A,2) /= size(c)) STOP "ERROR: Wrong dimensions for input arguments"

n = size(A,2)
x_p = 1./n
f_p = HUGE(f_p)
f = Potential(c, x_p)
delta = 1./8
eps = 1.E-4

IF( ANY(GEMV(A,x_p) >= eps) ) STOP "ERROR: Center of Simplex doesn't lie in Null space of the Input Matrix"

DO WHILE( f - f_p <= -delta)

	x = Optimize(x_p)
	f_p = f
	f =  Potential(c,x)	
	x_p = x

END DO

IF( DOT(c,x) >= eps )	STOP "ERROR: Optimum value is not zero."

x_can = x_p

CONTAINS
 
FUNCTION Optimize(x_p) RESULT(x)

REAL, INTENT(IN) :: x_p(:)
REAL :: x(size(x_p))
REAL :: e(size(x_p)), D(size(x_p),size(x_p)), AD(size(A,1),size(x_p)), B(size(A,1)+1,size(x_p)), I(size(x_p),size(x_p)), &
		 c_p(size(c)), c_unit(size(c)), x0(size(x_p)), alpha, r


e = 1.
D = Diag(x_p)
AD = GEMM(A,D)
B = AD .VAUG. e
I = Diag(e)

c_p = GEMV(GEMM((I - GEMM(TRANSPOSE(B),GEMM(SYMINV(GEMM(B, TRANSPOSE(B))),B))),D),c)

c_unit = c_p/ENORM(c_p)

x0 = 1./n
alpha = 1./4
r = 1./sqrt(real(n)*(real(n)-1))

x = x0 - alpha*r*c_unit

x = GEMV(D,x)/SUM(GEMV(D,x))

END FUNCTION Optimize

END FUNCTION Karmarkar
