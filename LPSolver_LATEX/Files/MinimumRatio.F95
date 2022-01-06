UREAL PURE FUNCTION minratio(n, cunit) RESULT(alpha)

INTEGER, INTENT(IN) :: n
UREAL, INTENT(IN) :: cunit(n)

UREAL, PARAMETER :: beta = ukind(1.Q-1)


alpha = (one-beta)/(n*maxval(cunit))

END FUNCTION minratio
