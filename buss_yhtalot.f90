MODULE buss_yhtalot
!Moduuli sisältää busselaattorin määrittelevät diff. yhtälöt Dx ja Dy.

CONTAINS

 FUNCTION Dx(x, y, t, a, b, alfa, o) RESULT (xtulos)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x, y, t, a, b, alfa, o
  REAL :: xtulos
  xtulos = a - (b + 1) * x + (x**2) * y + alfa * COS(o * t)
 END FUNCTION Dx

 FUNCTION Dy(x, y, t, a, b, alfa, o) RESULT (ytulos)
  IMPLICIT NONE
  REAL, INTENT(IN) :: x, y, t, a, b, alfa, o
  REAL :: ytulos
  ytulos = b * x - (x**2) * y
 END FUNCTION Dy

END MODULE buss_yhtalot
