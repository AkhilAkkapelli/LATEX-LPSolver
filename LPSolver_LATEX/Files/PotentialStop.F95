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
