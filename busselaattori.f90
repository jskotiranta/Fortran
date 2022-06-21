PROGRAM busselaattori

USE runge_kutta
USE buss_yhtalot

 IMPLICIT NONE
 REAL :: a, b, alfa, o, x1, y1, t1, xmin, xmax, ymin, ymax, tmin, tmax
 INTEGER :: tulostus
 INTEGER, PARAMETER :: pistlkm = 100000
 REAL, DIMENSION(pistlkm) :: xtau, ytau, ttau

 !Pyydet‰‰n parametrien arvot
 WRITE(*,*) 'Anna A'
 READ(*,*) a
 WRITE(*,*) 'Anna B'
 READ(*,*) b
 WRITE(*,*) 'Anna alfa'
 READ(*,*) alfa
 WRITE(*,*) 'Anna omega'
 READ(*,*) o
 WRITE(*,*) 'Anna x1 ja y1'
 READ(*,*) x1, y1
 WRITE(*,*) 'Anna tmin ja tmax.'
 READ(*,*) tmin, tmax
 t1 = tmin

 !Valitaan miss‰ tasossa halutaan tulostus n‰kyv‰ksi
 WRITE(*,*) 'Miss‰ tasossa k‰yr‰ piirret‰‰n?'
 WRITE(*,*) '1: (x,y)-taso.'
 WRITE(*,*) '2: (t,x)-taso.'
 WRITE(*,*) '3: (t,y)-taso.'
 WRITE(*,*) '4: Kaikki rinnakkain.'
 READ(*,*) tulostus
 IF (tulostus > 4) THEN
  WRITE(*,*) 'Virheellinen tulostusk‰sky.'
 END IF

 !Kutsutaan Runge-Kutan aliohjelmaa:
 CALL r_k4(Dx, Dy, a, b, alfa, o, x1, y1, t1, tmax, pistlkm, xtau,&
           & ytau, ttau)

 !Akselien minimi ja maksimiarvot, jotta kuvaaja n‰kyy j‰rkev‰sti.
 !Lis‰‰m‰ll‰ xmin:iin ja ymin:iin pieni termi v‰ltet‰‰n tilanne,
 !jossa minimi = maksimi eik‰ kuvaajaa piirry.
 xmin = MINVAL(xtau) - ((MAXVAL(xtau) - MINVAL(xtau))/10) - 0.0001
 xmax = MAXVAL(xtau) + ((MAXVAL(xtau) - MINVAL(xtau))/10)
 ymin = MINVAL(ytau) - ((MAXVAL(ytau) - MINVAL(ytau))/10) - 0.0001
 ymax = MAXVAL(ytau) + ((MAXVAL(ytau) - MINVAL(ytau))/10)
  

 !Tulostus pgplotilla:
 SELECT CASE(tulostus)
 CASE(1)
  !Tulostus (x,y)-tasossa:
  CALL pgopen('?')
  CALL pgenv(xmin, xmax, ymin, ymax, 0, 0)
  CALL pglabel('x','y','Busselaattori (x,y)-tasossa')
  CALL pgsci(8)
  CALL pgline(pistlkm, xtau, ytau)
  CALL pgsci(2)
  CALL pgpt1(xtau(1), ytau(1), 16)
  CALL pgend 
 CASE(2)
  !Tulostus (t,x)-tasossa:
  CALL pgopen('?')
  CALL pgenv(tmin, tmax, xmin, xmax, 0,0)
  CALL pglabel('t','x','Busselaattori (t,x)-tasossa')
  CALL pgsci(3)
  CALL pgline(pistlkm, ttau, xtau)
  CALL pgend
 CASE(3)
  !Tulostus (t,y)-tasossa:
  CALL pgopen('?')
  CALL pgenv(tmin, tmax, ymin, ymax, 0, 0)
  CALL pglabel('t', 'y' ,'Busselaattori (t,y)-tasossa')
  CALL pgsci(6)
  CALL pgline(pistlkm, ttau, ytau )
  CALL pgend
 CASE(4)
  CALL pgopen('?')
  CALL pgsubp(2,2)
  CALL pgsch(2.0)
  CALL pgenv(xmin, xmax, ymin, ymax, 0,0)
  CALL pglabel('x', 'y', ' ')
  CALL pgsci(8)
  CALL pgline(pistlkm, xtau, ytau)
  CALL pgsci(2)
  CALL pgpt1(xtau(1), ytau(1), 16)
  CALL pgsci(1)

  CALL pgenv(tmin, tmax, xmin, xmax, 0,0)
  CALL pglabel('t', 'x', ' ')
  CALL pgsci(3)
  CALL pgline(pistlkm, ttau, xtau)
  CALL pgsci(1)

  CALL pgenv(tmin, tmax, ymin, ymax, 0,0)
  CALL pglabel('t', 'y', ' ')
  CALL pgsci(6)
  CALL pgline(pistlkm, ttau, ytau)
  CALL pgsci(1)
  CALL pgend
 END SELECT

END PROGRAM busselaattori
