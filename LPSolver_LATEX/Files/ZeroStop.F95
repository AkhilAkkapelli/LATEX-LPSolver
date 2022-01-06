LOGICAL PURE FUNCTION zerostop(n, x,xp, iter, c) RESULT(stp)

INTEGER, INTENT(IN) :: iter, n
UREAL, INTENT(IN) :: x(n), xp(n), c(n)

INTEGER, PARAMETER :: iterlimit = 10000


stp = .FALSE.

IF(iter >= iterlimit .OR. x(n-1) < real(1.Q-100,realKind) .OR. x(n-1) > xp(n-1)) stp = .TRUE.

END FUNCTION zerostop
