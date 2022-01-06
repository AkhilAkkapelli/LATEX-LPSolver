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
