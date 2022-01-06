UREAL PURE FUNCTION potentialratio(n, cunit) RESULT(alpha)

INTEGER, INTENT(IN) :: n
UREAL, INTENT(IN) :: cunit(n)


alpha = 1/(4*NSQRT(ukind(n)*(ukind(n)-1)))

END FUNCTION potentialratio
