PROGRAM laneemden
IMPLICIT NONE

REAL :: n,ksi1,psi1,theta1,thetamin,thetamax,ksimin,ksimax,rhomin,rhomax,pmin,pmax,tmin,tmax,mmin,mmax,rho0,p0,t0,K,gamma,massa,sade,
INTEGER, PARAMETER :: pistlkm = 10000
REAL, PARAMETER :: Pi = 3.1415, G = 6.67259E-11
REAL, DIMENSION(pistlkm) :: ksitau,psitau,thetatau,rhotau,ptau,ttau,mtau
WRITE(*,*) 'Syötä n:n arvo'
READ(*,*) n
WRITE(*,*) 'Syota tiheys keskustassa'
READ(*,*) rho0
WRITE(*,*) 'Syötä paine keskustassa'
READ(*,*) p0
WRITE(*,*) 'Syötä tähden massa'
READ(*,*) massa
WRITE(*,*) 'Syötä lämpötila keskustassa'
READ(*,*) t0
ksi1 = 0.0001; psi1 = 0.0001; theta1 = 1.0

CALL r_k4(n,ksi1,psi1,theta1,pistlkm,ksitau,psitau,thetatau)

!Suureiden yhteydet thetaan
gamma = (1 / (n+0.001)) + 1
!K = ptau(1) / rhotau(1)**gamma
!sade = ((n + 1)*K/(4*Pi*G)) * rho0**((1-n)/(2*(n+0.001))) * ksitau(10000)
!p0 = G * massa**2 / sade**4 * (4 * Pi * (n+1) * psitau(10000)**2)**(-1)
rhotau = rho0 * thetatau**n
ptau = p0 * thetatau**(n + 1)
ttau = t0 * thetatau
!mtau!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 ksimin = 0.0001
 ksimax = Pi
 thetamin = MINVAL(thetatau) - ((MAXVAL(thetatau) - MINVAL(thetatau))/10) - 0.0001
 thetamax = MAXVAL(thetatau) + ((MAXVAL(thetatau) - MINVAL(thetatau))/10)
 rhomin = MINVAL(rhotau) - ((MAXVAL(rhotau) - MINVAL(rhotau))/10) - 0.001
 rhomax = MAXVAL(rhotau) + ((MAXVAL(rhotau) - MINVAL(rhotau))/10)
 pmin = MINVAL(ptau) - ((MAXVAL(ptau) - MINVAL(ptau))/10) - 0.0001
 pmax = MAXVAL(ptau) + ((MAXVAL(ptau) - MINVAL(ptau))/10)
 tmin = MINVAL(ttau) - ((MAXVAL(ttau) - MINVAL(ttau))/10) - 0.0001
 tmax = MAXVAL(ttau) + ((MAXVAL(ttau) - MINVAL(ttau))/10)
 mmin = MINVAL(mtau) - ((MAXVAL(mtau) - MINVAL(mtau))/10) - 0.0001
 mmax = MAXVAL(mtau) + ((MAXVAL(mtau) - MINVAL(mtau))/10)

  !Tulostus pgplottiin
  CALL pgopen('?')
  CALL pgsubp(2,2)
  CALL pgsch(2.0)

  CALL pgenv(ksimin, ksimax, thetamin, thetamax, 0, 0)
  CALL pglabel('\gc','\gh','Lane-Emden')
  CALL pgsci(8)
  CALL pgline(pistlkm, ksitau, thetatau)
  CALL pgsci(2)
  CALL pgpt1(ksitau(1), thetatau(1), 16)
  CALL pgsci(1)
  
  CALL pgenv(ksimin, ksimax, rhomin, rhomax, 0,0)
  CALL pglabel('\gc', '\gr', ' ')
  CALL pgsci(4)
  CALL pgline(pistlkm, ksitau, rhotau)
  CALL pgsci(2)
  CALL pgpt1(ksitau(1), rhotau(1), 16)
  CALL pgsci(1)

  CALL pgenv(ksimin, ksimax, pmin, pmax, 0,0)
  CALL pglabel('\gc', 'P', ' ')
  CALL pgsci(3)
  CALL pgline(pistlkm, ksitau, ptau)
  CALL pgsci(1)

  CALL pgenv(ksimin, ksimax, tmin, tmax, 0,0)
  CALL pglabel('\gc', 'T', ' ')
  CALL pgsci(6)
  CALL pgline(pistlkm, ksitau, ttau)
  CALL pgsci(1)
  CALL pgend

CONTAINS

 REAL FUNCTION Dtheta(theta, psi, ksi, n) RESULT (Dthetatulos)
  IMPLICIT NONE
  REAL :: ksi, psi, theta, n
  Dthetatulos = psi
 END FUNCTION Dtheta

 REAL FUNCTION Dpsi(theta, psi, ksi, n) RESULT (Dpsitulos)
  IMPLICIT NONE
  REAL :: ksi, psi, theta, n
  Dpsitulos = -2.0 * psi / ksi - theta**n
 END FUNCTION Dpsi

!Runge-Kutta iskee jälleen
SUBROUTINE r_k4(n,ksi1,psi1,theta1,pistlkm,ksitau,psitau,thetatau)
  IMPLICIT NONE
  REAL, INTENT(IN) :: n, ksi1, psi1, theta1
  INTEGER, INTENT(IN) :: pistlkm
  REAL, DIMENSION(pistlkm), INTENT(OUT) :: ksitau, psitau, thetatau
  REAL :: askel, kx1, kx2, kx3, kx4, ky1, ky2, ky3, ky4
  INTEGER :: i

  !Taulukoiden alustukset ja algoritmin askelkoko
  ksitau(1) = ksi1 ; psitau(1) = psi1; thetatau(1) = theta1
  askel = Pi / 10000

  !Neljännen kertaluvun Runge-Kutta algoritmi
  DO i = 1, (pistlkm - 1)
  ksitau(i + 1) = ksitau(i) + askel

   !Ensimmäiset gradientit saadaan i:nnestä "alkuarvosta":
   kx1 = askel * Dtheta(thetatau(i), psitau(i), ksitau(i), n)
   ky1 = askel * Dpsi(thetatau(i), psitau(i), ksitau(i), n)

   kx2 = askel * Dtheta(thetatau(i) + kx1/2.0, psitau(i) + ky1/2.0, ksitau(i) + askel/2.0, n)
   ky2 = askel * Dpsi(thetatau(i) + kx1/2.0, psitau(i) + ky1/2.0, ksitau(i) + askel/2.0, n)

   kx3 = askel * Dtheta(thetatau(i) + kx2/2.0, psitau(i) + ky2/2.0, ksitau(i) + askel/2.0, n)
   ky3 = askel * Dpsi(thetatau(i) + kx2/2.0, psitau(i) + ky2/2.0, ksitau(i) + askel/2.0, n)

   kx4 = askel * Dtheta(thetatau(i) + kx3, psitau(i) + ky3, ksitau(i) + askel, n)
   ky4 = askel * Dpsi(thetatau(i) + kx3, psitau(i) + ky3, ksitau(i) + askel, n)

   !Lopuksi gradienteistä saadaan lasketuksi (i + 1):nnet "alkuarvot":
   thetatau(i + 1) = (thetatau(i) + (kx1 + 2.0 * kx2 + 2.0 * kx3 + kx4) / 6.0)
   psitau(i + 1) = (psitau(i) + (ky1 + 2.0 * ky2 + 2.0 * ky3 + ky4) / 6.0)
  END DO

 END SUBROUTINE r_k4




END PROGRAM laneemden
