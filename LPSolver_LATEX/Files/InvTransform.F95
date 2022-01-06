FUNCTION InvTransform(x_art,a0) RESULT(x)

REAL, INTENT(IN) ::  x_art(:)
REAL, OPTIONAL :: a0(:)

REAL :: x(size(x_art)-1)

REAL :: x0(size(x_art)-1)


IF(.NOT. PRESENT(a)) THEN
	x0 = 1.
ELSE 
	x0 = a0
END IF

x = (x0*x_art)/x_art(size(x_art))

END FUNCTION InvTransform
