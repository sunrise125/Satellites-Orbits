      subroutine uv (r2, taux,sigma2,alfa,ff,gg )
  
      implicit none
      REAL(8), PARAMETER :: DPI = 3.141592653589793238462643D0    ! Pi
      Integer Sg
      REAL(8):: psi,psia,psir,psi0,trr,b,det,alfa,taux
      REAL(8):: S3,S2,S1,S0,fs,fs1,ff,gg,r2,sigma2    

     
     det = taux
     psia = det/r2  !  Initial value of Psi 
     
     IF (det >= 0D0) THEN 
      sg = 1
     ELSE 
      sg = -1
     END IF                     ! sg = sign for negative times
                                !If multiple REVOLUTIONS in Elliptical Orbits 
     IF (alfa < 0D0) THEN 
      psir = 2D0 * DPI / SQRT(-alfa)
      trr  = 2D0 *DPI / ((-alfa) * SQRT(-alfa))
      det = MODULO(det, trr)
      psia = MODULO(psia, psir)
          IF (ABS(det) >= trr / 2D0) THEN 
             det = (ABS(det) - trr) * sg
             psia = (ABS(psia) - psir) * sg
          END IF
     END IF
    
     DO
      b = alfa * psia**2D0
      S3 = psia**3D0/6D0*((((((((b/342D0+1D0)*b/272D0+1D0)*b/210D0+1D0)*b/156D0+1D0)  &
            * b/110D0+1D0)*b/72D0+1D0)*b/42D0+1D0)*b/20D0+1D0)
      S2 = psia**2D0/2D0*((((((((b/306D0+1D0)*b/240D0+1D0)*b/182D0+1D0)*b/132D0+1D0)  &
            * b/90D0+1D0)*b/56D0+1D0)*b/30D0+1D0)*b/12D0+1D0)
      S1 = psia+alfa*S3
      S0 = 1D0+alfa*S2

      fs = r2 * S1 + sigma2 * S2 + S3 - det     ! Time t 
      fs1= r2 * S0 + sigma2 * S1 + S2           ! r at time t 
      psi = psia - fs / fs1 
                                          ! NEWTON Iteration 
      psi0 = ABS(psi-psia)
      psia = psi
      
      IF ( psi0 < 0.000000001 ) THEN
       GOTO 10
      ELSE
       GOTO 5
      END IF
5      continue
     END DO 
10    ff = 1D0 - S2 / r2
      gg = det - S3                       ! ff , gg parameter  
   
     END 


   subroutine pveq (x,y,z,vx,vy,vz)
   implicit NONE
   INTEGER I, NR    !, openstatus
   REAL(8), PARAMETER :: DPI = 3.141592653589793238462643D0    ! Pi
   REAL(8), PARAMETER :: AS2R = 4.848136811095359935899141D-6  ! Arcseconds to radians
   REAL(8), PARAMETER :: DR2D = 57.29577951308232087679815D0   ! Radians to degrees
   REAL(8), PARAMETER :: D2PI = 6.283185307179586476925287D0   ! 2Pi
   REAL(8), PARAMETER :: DD2R = 1.745329251994329576923691D-2  ! Degrees to radians
   REAL(8), PARAMETER :: DR2AS = 206264.8062470963551564734D0  ! Radians to arc seconds
   REAL(8), PARAMETER :: DR2S = 13750.98708313975701043156D0   ! Radians to seconds
   REAL(8), PARAMETER :: MU = 398600.4418D0
   REAL(8), PARAMETER :: RT = 6378.137D0
   REAL(8), PARAMETER :: kgauss = 1.239446678D-3
   REAL(8), PARAMETER :: uv0 = 7.905364437
   character*20 FileName  
   REAL (8):: a,e,inc,capom,omega,capm
   REAL (8):: x,y,z,vx,vy,vz, r,V
   REAL (8):: n1,n2,n3, N,  wx,wy,wz, ifact, p,q,  sincapom,coscapom
   REAL (8):: fx,fy,fz,  gx,gy,gz,  ex,ey,ez,  h,k,  e2,  XX,YY
   REAL (8):: b, FF,sinF,cosF,  lambda,  zita, sinzita,coszita,  EE,VV
   REAL (8):: L,sinL,cosL
!
      print*,""
      print*,"  --- Given Position Vector [km] & Velocity Vector [km/s] ---  "
      write(*, *), x,y,z
      write(*, *), vx,vy,vz
! ----------------------------------------------------
! Computing semimajor axis
      r= sqrt(x*x + y*y + z*z); V= sqrt(vx*vx + vy*vy + vz*vz) 
      a= 1d0/(2d0/r - V*V/MU)      

! Unit vector w
      n1= y*vz - z*vy;    n2= z*Vx - x*vz;   n3= x*vy - y*vx
      N= sqrt(n1*n1 + n2*n2 + n3*n3)
      wx= n1/N;  wy= n2/N;  wz= n3/N

! Motion factor (ifact)
      IF (wz<0) then
         ifact= -1d0
          else
         ifact= 1d0
      END IF

! Equinoctial Elements (p,q) -> Node + Incl
      p= wx/(1d0+ifact*wz);  q= -wy/(1d0+ifact*wz)
         
      sincapom = p/sqrt(p*p+q*q);   coscapom = q/sqrt(p*p+q*q)
      capom= atan2(sincapom,coscapom)
       if (capom<0d0) then
          capom=capom+D2PI
       endif

      inc= DPI * (1d0 - ifact)/2d0 + 2d0*ifact*atan(sqrt(p*p+q*q))


! Unit vectors (f,g)
       fx= (1d0-p*p+q*q)/(1d0+p*p+q*q);  fy= (2d0*p*q)/(1d0+p*p+q*q);  fz= (-2d0*ifact*p)/(1d0+p*p+q*q)        
       gx= (2d0*ifact*p*q)/(1d0+p*p+q*q);  gy= (1d0+p*p-q*q)/(1d0+p*p+q*q);  gz= (2d0*q)/(1d0+p*p+q*q)

! Unit vector eccentricity (e)
       ex= -x/r + (vy*wz - vz*wy)*N/MU
       ey= -y/r + (vz*wx - vx*wz)*N/MU
       ez= -z/r + (vx*wy - vy*wx)*N/MU     
            e= sqrt(ex*ex + ey*ey + ez*ez)

! Equinoctial elements (h,k)
       h= ex*gx + ey*gy + ez*gz
       k= ex*fx + ey*fy + ez*fz
           e2= sqrt(h*h + k*k)

! Satellite coordinates (X,Y)
       XX= x*fx + y*fy + z*fz
       YY= x*gx + y*gy + z*gz

! Auxiliar parameter (b)
       b= 1d0/(1d0+sqrt(1d0-h*h-k*k))

! Eccentric longitude (F)
       sinF= h + ((1d0-h*h*b)*YY - h*k*b*XX)/ (a*sqrt(1d0-h*h-k*k))
       cosF= k + ((1d0-k*k*b)*XX - h*k*b*YY)/ (a*sqrt(1d0-h*h-k*k))
       FF= atan2(sinF,cosF)
         if (FF<0d0) then
          FF=FF+D2PI
         endif

! Mean longitude (lambda)
       lambda= FF +h*cosF - k*sinF

! Auxiliar parameter (zita)
       sinzita= h/sqrt(h*h + k*k);  coszita= k/sqrt(h*h + k*k)
        zita= atan2(sinzita,coszita)
         if (zita<0d0) then
          zita=zita+D2PI
         endif

! Argument of perigee (omega)
       omega= zita - ifact*capom
         if (omega<0d0) then
          omega=omega+D2PI
         endif

! Mean anomaly
       capm= lambda - zita
        if (capm<0d0) then
          capm=capm+D2PI
         endif

! Eccentric anomaly (E)
       EE= FF - zita
        if (EE<0d0) then
          EE=EE+D2PI
         endif       

! True longitude (L)
       sinL= ((1d0-k*k*b)*sinF+h*k*b*cosF-h)/(1-h*sinF-k*cosF)
       cosL= ((1d0-h*h*b)*cosF+h*k*b*sinF-k)/(1-h*sinF-k*cosF)
        L= atan2(sinL,cosL)
         if (L<0d0) then
          L=L+D2PI
         endif

! True anomaly (v)
       VV= L - zita
        if (VV<0d0) then
          VV=VV+D2PI
         endif  
             
! 

      print*,""
      print*," ---------- Output Satellite Orbital Elements ---------------------------"
      print*,'     Semi-major axis= ', a, " km"
      print*,'   1-   Eccentricity= ', e, " dimensionless" 
      print*,'   2-   Eccentricity= ', e2, " dimensionless"
      print*,"         Inclination= ", inc*DR2D, " degs"
      print*,"         Node (RAAN)= ", capom*DR2D, " degs"
      print*," Argument of Perigee= ", omega*DR2D, " degs"
      print*,"        Mean anomaly= ", capm*DR2D, " degs"
      print*,"        True anomaly= ", VV*DR2D, " degs"
      print*," ------------------------------------------------------------------------"
      print*,""

      end
 






