MODULE runge_kutta
!Moduuli ratkaisee kahden yhtälön yhtälöparin käyttäen neljännen kertaluvun
!Runge-Kutta menetelmää.

USE buss_yhtalot

CONTAINS

 SUBROUTINE r_k4(Dx,Dy,a,b,alfa,o,x1,y1,t1,tmax,pistlkm,xtau,ytau,ttau)
  IMPLICIT NONE
  REAL, EXTERNAL :: Dx, Dy
  REAL, INTENT(IN) :: a, b, alfa, o, x1, y1, t1, tmax
  INTEGER, INTENT(IN) :: pistlkm
  REAL, DIMENSION(pistlkm), INTENT(OUT) :: xtau, ytau, ttau
  REAL :: askel, kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4
  INTEGER :: i

  !Taulukoiden alustukset ja algoritmin askelkoko
  xtau(1) = x1 ; ytau(1) = y1; ttau(1) = t1
  askel = ((tmax - t1) / (REAL(pistlkm)))

  !Neljännen kertaluvun Runge-Kutta algoritmi
  DO i = 1, (pistlkm - 1)
   ttau(i + 1) = ttau(i) + askel

   !Ensimmäiset gradientit saadaan i:nnestä "alkuarvosta":
   kx1 = askel * Dx(xtau(i), ytau(i), ttau(i), a, b, alfa, o)
   ky1 = askel * Dy(xtau(i), ytau(i), ttau(i), a, b, alfa, o)

   !Toiset gradientit saadaan ensimmäisistä:
   kx2 = askel * Dx(xtau(i) + kx1/2.0, ytau(i) + ky1/2.0, ttau(i) + askel/2.0&
         &, a, b, alfa, o)
   ky2 = askel * Dy(xtau(i) + kx1/2.0, ytau(i) + ky1/2.0, ttau(i) + askel/2.0&
         &, a, b, alfa, o)

   !Kolmannet gradientit:
   kx3 = askel * Dx(xtau(i) + kx2/2.0, ytau(i) + ky2/2.0, ttau(i) + askel/2.0&
         &, a, b, alfa, o)
   ky3 = askel * Dy(xtau(i) + kx2/2.0, ytau(i) + ky2/2.0, ttau(i) + askel/2.0&
         &, a, b, alfa, o)

   !Neljännet gradientit:
   kx4 = askel * Dx(xtau(i) + kx3, ytau(i) + ky3, ttau(i) + askel, a, b, alfa,&
         & o)
   ky4 = askel * Dy(xtau(i) + kx3, ytau(i) + ky3, ttau(i) + askel, a, b, alfa,&
         & o)

   !Lopuksi gradienteistä saadaan lasketuksi (i + 1):nnet "alkuarvot":
   xtau(i + 1) = (xtau(i) + (kx1 + 2 * kx2 + 2 * kx3 + kx4) / 6.0)
   ytau(i + 1) = (ytau(i) + (ky1 + 2 * ky2 + 2 * ky3 + ky4) / 6.0)

  END DO

 END SUBROUTINE r_k4

END MODULE runge_kutta
