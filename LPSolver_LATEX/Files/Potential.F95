FUNCTION Potential(c,x) RESULT(f)

REAL, INTENT(IN) :: c(:) , x(:)
REAL :: f


IF(DOT(c,x) < 0) STOP 'Optimum value went negative: ', DOT(c,x)

f = SUM(LOG(DOT(c,x)/x))

END FUNCTION Potential
