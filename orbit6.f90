! -----------------------------------------------------------------------------------
!  PROGRAM orbit6.f90  -> Satellite Orbital Elements through 3 optical observations
!                         refined through Universal Variables - Canonical Units used        
!
!  Code in Fortran 90 + SOFA library written in Fortran 77 (222.for) + iers.for routine
!          (working folder should contain dbase files: iers.for, finals-iers.txt, obs.dbase like 447.dat) 
  
!          Observation data taken from p.447-449 of the book:
!          "Fundamentals of Astrodynamics and Applications" - 3rd ed. (2007) 
!          by David A. Vallado (ISBN 978-1-881883-14-2)
!
!  Compilation:
!          gfortran w.f90 1818.for iers.for twosubs.f90 -o w.exe  (Windows)
! 
!                        
!  Starting date: Apr.21,2013
!  Work in progress: Apr.23  (ended) 
! ------------------------------------------------------------------------------------
!
!--------------------------------------
!Date: 2012, 8, 20
!Location: lambda = -110 degs
!             phi =  +40 degs
!          height = 2000 m   
!-----------------------------
!Delta AT  to  internal routine
!DUT1      to  IERS.FOR routine
!Xp        to  IERS.FOR routine
!Yp        to  IERS.FOR routine 
!LOD       to  IERS.FOR routine
!----------------------
!  3 sat's observations (dbase)
!  ---------------------------------------
!  EXAMPLE 
!# 	Time (UTC) - August 20, 2012
!                    RA 	     RA(degs)       Decl        Decl(degs)
!	hh:mm:ss   	hh mm ss.ss 	     o 	      o  '  " 	       o
!1 	11:40:28   	00  03  45.58 	 0.939913 	+18 40 03.78 	+18.667717
!2 	23:48:28   	03  00  06.18 	45.025748 	+35 39 53.07 	+35.664741
!3 	23:52:28   	04  31  32.80 	67.886655 	+36 59 47.70 	+36.996553
! -------------------------------------------------------------------------------------
  PROGRAM orbitsei

    implicit NONE
   REAL(8), PARAMETER :: DPI = 3.141592653589793238462643D0    ! Pi
   REAL(8), PARAMETER :: AS2R = 4.848136811095359935899141D-6  ! Arcseconds to radians
   REAL(8), PARAMETER :: DR2D = 57.29577951308232087679815D0   ! Radians to degrees
   REAL(8), PARAMETER :: D2PI = 6.283185307179586476925287D0   ! 2Pi
   REAL(8), PARAMETER :: DD2R = 1.745329251994329576923691D-2  ! Degrees to radians
   REAL(8), PARAMETER :: DR2AS = 206264.8062470963551564734D0  ! Radians to arc seconds
   REAL(8), PARAMETER :: DR2S = 13750.98708313975701043156D0   ! Radians to seconds
   REAL(8), PARAMETER :: MU = 398600.4418D0                    ! Earth gravitational parameter: Km^3/s^2
   REAL(8), PARAMETER :: RT = 6378.137D0                       ! Conversion position factor: from canonical units to km
   REAL(8), PARAMETER :: MUc = 1D0                             ! Canonical gravit. parameter
   REAL(8), PARAMETER :: kgauss = 0.0012394466678D0            ! Earth central motion gravit. param. (times in seconds)
   REAL(8), PARAMETER :: uv0 = 7.905364437D0                   ! Conversion velocity factor: from canonical units to km/s
   
   INTEGER  IY,MO,ID,H1(3),MI1(3),J,I,SGN,H(3),M(3),FINAL,NR,kount,iora(4),iora2(4),iora3(4)
   INTEGER  openstatus,aux1,aux2
   REAL(8):: longit,latit,LON,LAT,height,SEC1(3),S(3),T(3),DJMJD0,LOD(3),UTC(3),TIME
   REAL(8):: UTC1(3), UTC2(3), UTC3(3),UT11(3), UT12(3), UT(3),TAI1(3),TAI2(3)
   REAL(8):: DATE(3),D(3),M0(3),S0(3),DEC(3),TUT(3),L(3,3),L0,LL(3,3),U,V,XYZ(3),RA(3),C,CTY,JD,EPS,GN
   REAL(8):: XP(3),YP(3),DUT1(3),DX0(3),DY0(3),Robs,L1(3),L2(3),L3(3),Ts,Tm,Np
   REAL(8):: iau_ANP,iau_GMST00,iau_S06,GST(3),LAST(3),DX06(3),DY06(3)
   REAL(8):: iau_ERA00, iau_EE00,iau_SP00,DDP00,DDE00
   REAL(8):: X,Y,XP0(3),YP0(3),TT1(3),TT2(3),TT(3),r(3,3),RECI(3,3)
   REAL(8):: RR1(3),RR2(3),RR3(3),r1,r2,r3,Vv,ei,ej,ek,ecc,dI,dJ,dK
   REAL(8):: ialpha,capom,omega,capm
! -------------------------------------------------------   
   REAL(8):: Vect1(3),Vect2(3),Vect3(3),mVe1,mVe2,mVe3,time12,theta12 ,sigma2 ,alfa
   REAL(8):: zz, zit, zzz, Cot, amin, Ag, B1, B3, aprec, psia, psi, psi0,r3c
   REAL(8):: Sf,Cf, Spri,Cpri, yy,xx, F0,F1, r1c,r2c, b, S1, S2, S3, S00
   REAL(8):: aorb,porb,eorb, V1orb,V2orb,V3orb,ff,gg,gg1,ff1,ff3,gg3,V2,rprim
   REAL(8):: V1x,V1y,V1z, V2x,V2y,V2z,x1,x2,x3,y3,z3,y1,y2,z1,z2,period,tp,taux
   REAL(8):: Vx2,Vy2,Vz2, dn1,dn3,fg
! ---------------------------------------------------   
   REAL(8):: dd1(3),dd2(3),dd3(3), nn1(3),nn2(3),nn3(3),RHO(3), Rhoaux(3), CC(3) 
   REAL(8):: tau1,tau3, ttt, aux, a, a1,a3, a1u,a3u, d1,d2, Lt,Mt,Kt, rsite2
   REAL(8):: traRC2IT1(3,3),traRC2IT2(3,3),traRC2IT3(3,3),xx1(3),xx2(3),xx3(3)
   REAL(8):: Yi,Yj,Yk,W,vx,vy,vz,px,py,pz,DetL,MM(3,3),Modec,rsec,a_prec
   REAL(8):: x1sto(3),x2sto(3),x3sto(3),xtemp(3)
   REAL(8):: ModD, ModN, prectum, semiax, ModS, Lscalar, ModB, ModV ,Nodo,Om    
   REAL(8):: incl,hi,hj,hk,Wi,Wj,WW,aW,Pw,n,Eccan,Man,argp,node,RPOM(3)
   REAL(8):: Ltn,Mtn,Ktn, f, fderiv, xl, deltax, uu, c1,c2,c3,Modh,Modnn
   REAL(8):: R12(3),R23(3),R31(3), R12N(3),R23N(3),R31N(3), R23S(3),R31S(3),R12S(3) 
   REAL(8):: Dg(3),Ng(3),Sg(3), Bg(3),Vg(3), hh(3),kk(3),nn(3), ec(3),vh(3)
   REAL(8):: Ss(3),X0(3),Y0(3),RC2I(3,3),ERA(3),RC2TI(3,3),RC2IT1(3,3),RC2IT2(3,3)
   REAL(8):: xy1(3),xy2(3),xy3(3),XYZkm(3),inerkm(3),INERTIAL(3),RC2IT3(3,3)

   character*20 FileName,tr 
   character*2 ch,FL,MA 
   character*1 SIGN        
! --------------------------------------------  
   REAL(8), dimension(:), allocatable :: STORE
! --------------------------------------------       
      WRITE(*,*)""
      FINAL = 58940       ! 2020/04/01   reported in MJD IERS final database (for update) 
      WRITE (*,*) " MJD currently final number, FINAL = ",FINAL," 2020/04/01 "
      WRITE(*,*)""
      WRITE (*,*) "  --- Starting value for Lagrange equation assumption: xl= 4.0 (canonical units) ---"
      WRITE (*,*) "   (to be increased in case the convergence deals to negative root)"

! --- Start program --------
      write(*, *)
      write(*, *)
      write(*, *) '=========================================================='
      write(*, *) 'Satellite Orbital Elements through 3 optical observations'
      write(*, *) '     (Gauss method, adapted from Vallado algorithms)'
      write(*, *) '=========================================================='
      write(*, *)
 

! ---------------- Start data file reading ------------- 
      WRITE (*, '(1X, A)', ADVANCE = "NO") "Enter name of data file: "
      READ *, FileName   
      OPEN (99, FILE = FileName, STATUS = "OLD", IOSTAT=OpenStatus)
      IF (OpenStatus > 0) STOP "   ***Cannot open the file ***"
      read(99, *), NR
      allocate(STORE(NR))
      read(99,*) STORE
! ---------------- End data file reading ------------- 
!      Checkup only
!      DO I=1,NR 
!      write(*, *), STORE(I)
!      END DO 

! -------- Assignment of stored values to own variables --------      
      longit= STORE(1)  ! Site eastern longitude (western, negativa)
      latit= STORE(2)   ! Site latitude
      height= STORE(3)  ! Site, meters above sea level
!  Division by 1 is a trick to transform numbers stored as real into integer    
      IY= STORE(4)/1  ! year
      MO= STORE(5)/1  ! month
      ID= STORE(6)/1  ! day

!     H1(I), MI1(I), SEC1(I) -> Time (hours, minutes, seconds) for 3 obs
      H1(1)= STORE(7);      MI1(1)= STORE(8);      SEC1(1)= STORE(9)
      H1(2)= STORE(16);     MI1(2)= STORE(17);     SEC1(2)= STORE(18)
      H1(3)= STORE(25);     MI1(3)= STORE(26);     SEC1(3)= STORE(27)  
        
!     H(I), M(I), S(I) -> RA (hours, minutes, seconds) for 3 obs
      H(1)= STORE(10);      M(1)= STORE(11);      S(1)= STORE(12)
      H(2)= STORE(19);      M(2)= STORE(20);      S(2)= STORE(21) 
      H(3)= STORE(28);      M(3)= STORE(29);      S(3)= STORE(30)       
      
!     D(I), M0(I), S0(I) -> Decl (degs, mins, arcsecs) for 3 obs
      D(1)= STORE(13);      M0(1)= STORE(14);      S0(1)= STORE(15)      
      D(2)= STORE(22);      M0(2)= STORE(23);      S0(2)= STORE(24)      
      D(3)= STORE(31);      M0(3)= STORE(32);      S0(3)= STORE(33)       
!
! ------------------------------------------------      
      print*,""
      print*,"    ----- DATE and LOCATION -----"
      WRITE(*,100) IY,MO,ID  
 100  FORMAT(" Date (YYYY:MM:DD)= ", i5,i3,i3)
      WRITE(*,101) longit,latit,height  
 101  FORMAT(" Longit(degs), Latit(degs), Height(m)= ", 2f11.6,f8.0)
 
      print*,""
      print*,"      ----- 1st OBSERVATION -----"
      WRITE(*,105) H1(1),MI1(1),SEC1(1)  
 105  FORMAT(" Time (HH:MM:SS.SS)= ", i3,i4,f7.2)
      WRITE(*,107) H(1),M(1),S(1) 
 107  FORMAT("   RA (HH:MM:SS.SS)= ", i3,i4,f7.2)
      aux1=STORE(13)/1; aux2=STORE(14)/1;          ! trick to print integers
      WRITE(*,108) aux1, aux2 , S0(1)  
 108  FORMAT("  Dec (dd:mm:ss.ss)= ", i3,i4,f7.2)
      print*,""

      print*,"      ----- 2nd OBSERVATION -----"
      WRITE(*,205) H1(2),MI1(2),SEC1(2)  
 205  FORMAT(" Time (HH:MM:SS.SS)= ", i3,i4,f7.2)
      WRITE(*,207) H(2),M(2),S(2) 
 207  FORMAT("   RA (HH:MM:SS.SS)= ", i3,i4,f7.2)    
      aux1=STORE(22)/1; aux2=STORE(23)/1;          ! trick to print integers
      WRITE(*,208) aux1, aux2 , S0(2)  
 208  FORMAT("  Dec (dd:mm:ss.ss)= ", i3,i4,f7.2)
      print*,""

      print*,"      ----- 3rd OBSERVATION -----"
      WRITE(*,305) H1(3),MI1(3),SEC1(3)  
 305  FORMAT(" Time (HH:MM:SS.SS)= ", i3,i4,f7.2)
      WRITE(*,307) H(3),M(3),S(3) 
 307  FORMAT("   RA (HH:MM:SS.SS)= ", i3,i4,f7.2)
      aux1=STORE(31)/1; aux2=STORE(32)/1;          ! trick to print integers
      WRITE(*,308) aux1, aux2 , S0(3)  
 308  FORMAT("  Dec (dd:mm:ss.ss)= ", i3,i4,f7.2)
      print*,""
!------------------------------------------------------------------------------------ 
!      
! Times in seconds 
      DO I = 1,3
        T(I) = H1(I) * 3600D0 + MI1(I) * 60D0 + SEC1(I)
      END DO       
! Right Ascensions (RA) in radians      
      DO I = 1,3
       RA(I) = ((H(I) + M(I)/60D0 + S(I)/3600D0) * 15D0 ) * DD2R     
      END DO   
     
! Declinations (Dec) in radians          
      DO I = 1,3
       IF( D(I) < 0D0 ) THEN
        DEC(I) = ((D(I) - M0(I)/60D0 - S0(I)/3600D0)) * DD2R
       ELSE IF ( D(I) >= 0D0 ) THEN
       DEC(I) = ((D(I) + M0(I)/60D0 + S0(I)/3600D0)) * DD2R
       END IF     
      END DO    


      LAT = latit * DD2R            ! latitude in rad.
      LON = longit * DD2R           ! longitude in rad.
      DJMJD0 = 2400000.5D0
      
      
!---------STEP 3 form matrix
      DO I = 1,3
         L1(I) = COS(RA(I)) * COS(DEC(I))                 
         L2(I) = SIN(RA(I)) * COS(DEC(I))
         L3(I) = SIN(DEC(I))
      END DO

      L(1,1) = L1(1);       L(2,1) = L2(1);      L(3,1) = L3(1)
      L(1,2) = L1(2);       L(2,2) = L2(2);      L(3,2) = L3(2)
      L(1,3) = L1(3);       L(2,3) = L2(3);      L(3,3) = L3(3) 


!----------- L = determinant of matrix    

      DetL= L(1,1)*L(2,2)*L(3,3) + L(2,1)*L(3,2)*L(1,3) + L(1,2)*L(2,3)*L(3,1) &
           -L(1,3)*L(2,2)*L(3,1)- L(2,3)*L(3,2)*L(1,1)- L(1,2)*L(2,1)*L(3,3)
    
!      
      LL(1,1)= (L(2,2)*L(3,3)-L(2,3)*L(3,2))/DetL;   
               LL(2,1)=(-L(2,1)*L(3,3)+L(2,3)*L(3,1))/DetL;   
                        LL(3,1)= (L(2,1)*L(3,2)-L(2,2)*L(3,1))/DetL;
      LL(1,2)=(-L(1,2)*L(3,3)+L(3,2)*L(1,3))/DetL;   
               LL(2,2)= (L(1,1)*L(3,3)-L(1,3)*L(3,1))/DetL;   
                        LL(3,2)=(-L(1,1)*L(3,2)+L(1,2)*L(3,1))/DetL;
      LL(1,3)= (L(1,2)*L(2,3)-L(1,3)*L(2,2))/DetL;   
               LL(2,3)=(-L(1,1)*L(2,3)+L(1,3)*L(2,1))/DetL;   
                        LL(3,3)= (L(1,1)*L(2,2)-L(1,2)*L(2,1))/DetL;

 
!-----------CALL UTC
  
     DO I = 1,3
      CALL iau_DTF2D ( "UTC",IY,MO,ID,H1(I),MI1(I),SEC1(I),UTC1(I),UTC2(I),J )  ! UTC1 + UTC2 = Julian Day
      IF ( J.NE.0 ) STOP                                                        ! UTC1 is the Julian Day number and
                                                                                ! UTC2 is the fraction of a day.   
      UTC(I) = UTC1(I) + UTC2(I) - DJMJD0                                       ! UTC in MJD (modified JD)  
      DATE(I) = UTC1(I) - DJMJD0                                                ! DATE is the MJD number
     
!----------Call the interpolation routine for per XP,YP,DUT1,dX,dY
       
      CALL iers_calc(UTC(I),DATE(I),FINAL,XP0(I),YP0(I),DUT1(I),LOD(I),DX0(I),DY0(I))
                               
      UTC2(I) = UTC2(I) + LOD(I) /1000D0/86400D0 
      TUT(I)  = UTC2(I) + DUT1(I)/86400D0        !  in fraction of day          
     
      CALL iau_SXP (AS2R, XP0, XP)               !  arcsec---> radians
      CALL iau_SXP (AS2R, YP0, YP)   

      CALL iau_SXP (( AS2R / 1000D0),DX0,DX06 )  !  mas---> radians 
      CALL iau_SXP (( AS2R / 1000D0),DY0,DY06 )
     
!---------- UTC -> TAI -> TT -> TCG.
      CALL iau_UTCTAI ( UTC1(I), UTC2(I), TAI1(I), TAI2(I), J )
      IF ( J.NE.0 ) STOP
      CALL iau_TAITT ( TAI1(I), TAI2(I), TT1(I), TT2(I), J )
      IF ( J.NE.0 ) STOP

      TT(I) = TT1(I) + TT2(I) - DJMJD0        ! Dynamical time in MJD
    
     END DO
      

!-- Transform geodetic to geocentric coordinates of the site. ECEF (reference ellipsoid  WGS84)

  
       CALL iau_GD2GC ( 1, LON, LAT, height, XYZ, J )   ! Geodetic to Geocentric  coordinate
       
       IF ( J.NE.0 )STOP
       U = SQRT ( XYZ(1)*XYZ(1) + XYZ(2)*XYZ(2) )    ! U = distance from Earth spin axis (km)
       V = XYZ(3)                                    ! V = distance north of equatorial plane (km)
 
       Robs = (SQRT ( U*U + V*V))/1000D0             ! Topocentic distance from Earth's center              
       XYZ = XYZ / 1000D0

!========================================================================================================

 ! ---------------------   Start Code 1  --------------------------
  

!   CIP  and  CIO,  IAU  2006/2000A.
       CALL   iau_XY06   (  DJMJD0,   TT(1),  X, Y  ) 
        S =  iau_S06  (  DJMJD0,   TT(1),  X,  Y )

!   Add  CIP  corrections.
       X  = X  + DX06(1) 
       Y  = Y  + DY06(1) 

!   GCRS  to  CIRS  matrix.
       CALL   iau_C2IXYS  (   X,  Y,  S, RC2I   ) 

!   Earth  rotation  angle.
       ERA(1)  =  iau_ERA00  (   DJMJD0+DATE(1),  TUT(1)   )

!   Form  celestial-terrestrial  matrix   (no polar motion  yet).
       CALL   iau_CR   ( RC2I,   RC2TI  ) 
       CALL   iau_RZ   ( ERA(1),  RC2TI   ) 

!   Polar  motion   matrix   (TIRS->ITRS,  IERS     2003).
       CALL   iau_POM00   (  XP(1),  YP(1),  iau_SP00(DJMJD0,TT(1)),  RPOM )

!   Form  celestial-terrestrial  matrix  (including polar  motion).
       CALL   iau_RXR   ( RPOM,   RC2TI,   RC2IT1  ) 

! ......... Call for Degrees to Hours, Minutes, Seconds + Fraction
       CALL   iau_A2TF   ( 9,  ERA(1), sign, iora )
! ---------------------------
!  New calc.
!  -----  Inertial Vector -------- 
       CALL  iau_TR  (RC2IT1 , traRC2IT1)
       CALL   iau_RXP  (traRC2IT1, XYZ, xx1)  ! Polar Motion Vector in km
       CALL   iau_SXP    (1D0, XYZ, XYZkm)  ! Fixed Pole Vector in km      

! -------- copy xx1 into x1sto, storage matrix -----    
       CALL iau_CP (xx1,x1sto)

! ---------------------     End Code 1   --------------------------

! 
!==============================================================================


! ---------------------   Start Code 2  --------------------------
!   CIP  and  CIO,  IAU  2006/2000A.
       CALL   iau_XY06   (  DJMJD0,   TT(2),  X, Y  ) 
        S =  iau_S06  (  DJMJD0,   TT(2),  X,  Y )

!   Add  CIP  corrections.
       X  = X  + DX06(2) 
       Y  = Y  + DY06(2) 

!   GCRS  to  CIRS  matrix.
       CALL   iau_C2IXYS  (   X,  Y,  S, RC2I   ) 
   
!   Earth  rotation  angle.
       ERA(2)  =  iau_ERA00  (   DJMJD0+DATE(2),  TUT(2)   )

!   Form  celestial-terrestrial  matrix   (no polar motion  yet).
       CALL   iau_CR   ( RC2I,   RC2TI  ) 
       CALL   iau_RZ   ( ERA(2),  RC2TI   ) 

!   Polar  motion   matrix   (TIRS->ITRS,  IERS   2003).
       CALL   iau_POM00 ( XP(2), YP(2), iau_SP00(DJMJD0,TT(2)), RPOM )

!   Form  celestial-terrestrial  matrix  (including polar  motion).
       CALL   iau_RXR   ( RPOM,   RC2TI,   RC2IT2  ) 

! ......... Call for Degrees to Hours, Minutes, Seconds + Fraction
       CALL   iau_A2TF   ( 9,  ERA(2), sign, iora )
      
! ---------------------------
!  New calc.
!  -----  Inertial Vector --------  
       CALL   iau_TR     (RC2IT2, traRC2IT2)  ! use transpose matrix   
       CALL   iau_RXP    (traRC2IT2, XYZ, xx2)            
       CALL   iau_SXP    (1D0, XYZ, XYZkm)  ! Fixed Pole Vector in km

! -------- copy xx2 into x2sto, storage matrix -----    
       CALL iau_CP (xx2,x2sto)


! ---------------------     End Code 2  --------------------------
!===============================================================================
! ---------------------   Start Code 3  --------------------------

!   CIP  and  CIO,  IAU  2006/2000A.
       CALL   iau_XY06   (  DJMJD0,   TT(3),  X, Y  ) 
        S =  iau_S06  (  DJMJD0,   TT(3),  X,  Y )

!   Add  CIP  corrections.
       X  = X  + DX06(3) 
       Y  = Y  + DY06(3) 

!   GCRS  to  CIRS  matrix.
       CALL   iau_C2IXYS  (   X,  Y,  S, RC2I   ) 

!   Earth  rotation  angle.
       ERA(3)  =  iau_ERA00  (   DJMJD0+DATE(3),  TUT(3)   )

!   Form  celestial-terrestrial  matrix   (no polar motion  yet).
       CALL   iau_CR   ( RC2I,   RC2TI  ) 
       CALL   iau_RZ   ( ERA(3),  RC2TI   ) 

!   Polar  motion   matrix   (TIRS->ITRS,  IERS     2003).
       CALL   iau_POM00   (  XP(3),  YP(3),  iau_SP00(DJMJD0,TT(3)),  RPOM )

!   Form  celestial-terrestrial  matrix  (including polar  motion).
       CALL   iau_RXR   ( RPOM,   RC2TI,   RC2IT3  ) 

! ......... Call for Degrees to Hours, Minutes, Seconds + Fraction
       CALL   iau_A2TF   ( 9,  ERA(3), sign, iora )
! ---------------------------
!  New calc.
!  -----  Inertial Vector -------- 
       CALL   iau_TR     (RC2IT3, traRC2IT3)  ! use transpose matrix    
       CALL   iau_RXP    (traRC2IT3, XYZ, xx3)  ! Polar Motion Vector in km
       CALL   iau_SXP    (1D0, XYZ, XYZkm )  ! Fixed Pole Vector in km

! -------- copy xx3 (3) into xx3sto(3), storage matrix -----    
       CALL iau_CP (xx3,x3sto)

! ---------------------     End Code 3  --------------------------


 
!    Start output
! -------------------
     WRITE(*,*)"==================================================================================="
     WRITE(*,*)"                 ---------  IERS EOP DATA INPUT -------                                 "
     WRITE(*,*)
     WRITE(*,499),TT(1),TT(2),TT(3)
 499 FORMAT(" TT(1)",f18.12,4x,"TT(2)",f18.12,4x"TT(3)",f18.12)    
     WRITE(*,*) 
     WRITE(*,500)DUT1(1),DUT1(2),DUT1(3)
 500 FORMAT(" DUT1(1)",f10.7," sec.",5x,"DUT1(2)",f10.7," sec.",5x,"DUT1(3)",f10.7," sec.")
     WRITE(*,*)
     WRITE(*,501)LOD(1),LOD(2),LOD(3)
 501 FORMAT(" LOD(1)",2x,f6.4,4x,"msec.",4x,"LOD(2)",2x,f6.4,4x,"msec.",4x,"LOD(3)"2x,f6.4,4x,"msec.")
     WRITE(*,*)   
     WRITE(*,502)XP0(1),XP0(2),XP0(3)  
 502 FORMAT(" XP(1)",2x,f9.6,x," arcsec.",x," XP(2)",2x,f9.6,x," arcsec.",x," XP(3)",2x,f9.6,x," arcsec.") 
     WRITE(*,*)  
     WRITE(*,503),YP0(1),YP0(2),YP0(3)
 503 FORMAT(" YP(1)",2x,f9.6,x," arcsec.",x," YP(2)",2x,f9.6,x," arcsec.",x," YP(3)",2x,f9.6,x," arcsec.") 
     WRITE(*,*)"==================================================================================="
     WRITE(*,*)       



!-------- Earth  rotation  angle.
     WRITE(*,*)""
     WRITE(*,*) " ERA1 ",ERA(1) * DR2D,"degs" 
     WRITE(*,*) " ERA2 ",ERA(2) * DR2D,"degs" 
     WRITE(*,*) " ERA3 ",ERA(3) * DR2D,"degs"            
     WRITE(*,*)

!------- ECI Matrix can be determined using the following equations:

      RECI(1,1) = x1sto(1);     RECI(1,2) = x2sto(1);     RECI(1,3) = x3sto(1)
      RECI(2,1) = x1sto(2);     RECI(2,2) = x2sto(2);     RECI(2,3) = x3sto(2)
      RECI(3,1) = x1sto(3);     RECI(3,2) = x2sto(3);     RECI(3,3) = x3sto(3)         
  
! -------- New approach (via Universal Variables) ---------
!
     print*,"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
! ---- RECI canonical ----
      RECI(1,1) = x1sto(1)/RT;     RECI(1,2) = x2sto(1)/RT;     RECI(1,3) = x3sto(1)/RT
      RECI(2,1) = x1sto(2)/RT;     RECI(2,2) = x2sto(2)/RT;     RECI(2,3) = x3sto(2)/RT
      RECI(3,1) = x1sto(3)/RT;     RECI(3,2) = x2sto(3)/RT;     RECI(3,3) = x3sto(3)/RT          
      WRITE(*,*)
      WRITE(*,*),"----- ECI Site Vectors MATRIX [Canonical Units] -----"
      WRITE(*,*)""
      WRITE(*,*),RECI(1,1),RECI(1,2),RECI(1,3)
      WRITE(*,*),RECI(2,1),RECI(2,2),RECI(2,3)
      WRITE(*,*),RECI(3,1),RECI(3,2),RECI(3,3)
      WRITE(*,*)""

! -------- Computing the Matrix MT = L^(-1) * RECI ----  
     CALL iau_RxR (LL,RECI,MM)
     
! ------- Dot product C = L2*RECL2, related to 2nd obs. ---------- 
     C= L(1,2)*RECI(1,2)+ L(2,2)*RECI(2,2)+ L(3,2)*RECI(3,2)
     print*,' C [canonical units]= ',C  
     print*,""
! -------- Computing delta-times in Canonical Units ----------
     tau1= kgauss*(T(1)-T(2))    
     tau3= kgauss*(T(3)-T(2))  
  
! -------- Computing time parameters (a1,a3) & (a1u,a3u) & output ---------- 
     aux= tau3-tau1
       a1= tau3/aux
       a3= -tau1/aux   
       a1u= tau3*(aux**2-tau3**2)/(6*aux)
       a3u= -tau1*(aux**2-tau1**2)/(6*aux)
     print*,'     (a1,a3)= ', a1,a3 
     print*,'   (a1u,a3u)= ', a1u,a3u 

! -------- Computing parameters (d1,d2) & output ----------
     d1=MM(2,1)*a1 - MM(2,2) + MM(2,3)*a3
     d2=MM(2,1)*a1u + MM(2,3)*a3u
     print*, '     (d1,d2)= ', d1,d2
     print*,""
!--------------------------------------------    

! Distance Site to Earth's center [km],  from site coords at time t2                 
     rsite2 = SQRT(RECI(1,2)**2+RECI(2,2)**2+RECI(3,2)**2)  
     Lt= d1**2 + 2*C*d1 + rsite2**2
     Mt= 2*(C*d2 + d1*d2)
     Kt= d2**2
! -------- Computing Lagrange 8th degree coefficients (Lt,Mt,Kt) & output ----------
!     x=r2
!                x^8 - Lt*x^6 - MU*Mt*x^3 - MU^2*Kt = 0
!     being:
!            Lt= d1^2 + 2*C*d1 + rsite2^2  
!            Mt= 2*(C*d2+d1*d2)
!            kt= d2^2 
!                and
!            MU= 1.0  (Canonical Units)
! ----------------------------------------------------------------------

     print*,""
     print*,'      MU = ',MUc, "Geocentric gravitational constant (=1, in Canonical Units)" 
     print*,"     Lagrange equation ->  x^8 - Lt*x^6 - MU*Mt*x^3 - MU^2*Kt = 0"
     print*,'     (Lt)= ', Lt
     print*,'     (mu*Mt)= ', MUc*Mt
     print*,'     (mu^2*Kt)= ', MUc**2*Kt
     print*,""
!
      kount=0
      xl= 4d0   ! Starting value     
      DO     
          f = xl**8 - Lt*xl**6 - MUc*Mt*xl**3 - MUc**2*Kt
          fderiv= 8*xl**7 - 6*Lt*xl**5 - 3*MUc*Mt*xl**2
         deltax= f/fderiv
          xl = xl - deltax
            IF (ABS(deltax).lt. 1D-6) then
             exit
            END IF 
        kount=kount+1    
      END DO
      
      print*,'     r2= ', xl, ' -> iterations used = ', kount
      print*,""

! ----- Computing value uu=mu/(r2)^3 and coefficients c1,c2,c3 ------- 
      r2 = xl    
      uu = MUc/r2**3D0      
      
      ff1 = 1d0 - uu*tau1*tau1/2d0              
      ff3 = 1d0 - uu*tau3*tau3/2d0
      gg1 = tau1*(1d0-uu*tau1*tau1/6d0)         
      gg3 = tau3*(1d0-uu*tau3*tau3/6d0)
      
      fg = ff1*gg3 - ff3*gg1
      
      c1 =  gg3/fg                      
      c2 = -1D0
      c3 =  -gg1/fg                     

      print*,' (f3,g3)= ', ff3,gg3
!      print*,""  
      print*,' (f1,g1)= ', ff1,gg1
      print*,""    

! Storing to vector CC(i)------------------------------------------
        CC(1)= -c1
        CC(2)= -c2
        CC(3)= -C3 
     print*,""
     print*,'     (u)= ', uu    
     print*,' (c1,c3)= ', c1,c3
     print*,'    (c2)= ', c2
     print*,""
! -------------------------------------------------------------------

     kount=-1; a=0d0
1234 continue
     aprec=a ! previous value of iteration
     kount= kount+1
! ----- Computing RHO(i), topocentric vectors ---------- i=1,2,3 ----
     CALL iau_RXP(MM,CC,Rhoaux)
     RHO(1)= Rhoaux(1)/c1
     RHO(2)= Rhoaux(2)/c2
     RHO(3)= Rhoaux(3)/c3    
!
     print*, "--------- Iteration nr.", kount," ----------"
     print*, '  (rho1)= ', RHO(1)*RT," km"
     print*, '  (rho2)= ', RHO(2)*RT," km"
     print*, '  (rho3)= ', RHO(3)*RT," km"
     print*, '  ' 

! ----------- Geocentric Vectors -------       
         RR1(1)= RHO(1)* L(1,1) + RECI (1,1)
         RR1(2)= RHO(1)* L(2,1) + RECI (2,1)
         RR1(3)= RHO(1)* L(3,1) + RECI (3,1)
!           
         RR2(1)= RHO(2)* L(1,2) + RECI (1,2)    
         RR2(2)= RHO(2)* L(2,2) + RECI (2,2)
         RR2(3)= RHO(2)* L(3,2) + RECI (3,2)        
!
         RR3(1)= RHO(3)* L(1,3) + RECI (1,3)    
         RR3(2)= RHO(3)* L(2,3) + RECI (2,3)
         RR3(3)= RHO(3)* L(3,3) + RECI (3,3) 

      x1 = RR1(1);  x2 = RR2(1);  x3 = RR3(1)
      y1 = RR1(2);  y2 = RR2(2);  y3 = RR3(2)
      z1 = RR1(3);  z2 = RR2(3);  z3 = RR3(3)      
                 
! ---- 3 geocentric vectors MODULUS (r1,r2,r3)     
      CALL iau_PM(RR1,r1)
      CALL iau_PM(RR2,r2)
      CALL iau_PM(RR3,r3)   
!                 
!      print*,""
!      print*,"----- GEOCENTRIC Vectors [RT] --------"
!      print*,""
!      print*,"r1",RR1(1),RR1(2),RR1(3)
!      print*,"r2",RR2(1),RR2(2),RR2(3)
!      print*,"r3",RR3(1),RR3(2),RR3(3)
!      print*,"" 
       print*,'   Modulus-r1= ', r1*RT, " km"
       print*,'   Modulus-r2= ', r2*RT, " km" 
       print*,'   Modulus-r2= ', r3*RT, " km"
       print*,"" 
            
!-------- 3 geocentric velocity vector ( V2x,V2y,V2z)
      dn1 = -ff3/fg;     dn3 = ff1/fg; 
!
      V2x = dn1 * x1 + dn3 * x3 
      V2y = dn1 * y1 + dn3 * y3
      V2z = dn1 * z1 + dn3 * z3
!--------  geocentric velocity modulus (V2)
      V2 = SQRT(V2x**2D0 + V2y**2D0 + V2z**2D0)

     print*,"   Velocity-V2x= ", V2x*uv0, "km/s"
     print*,"   Velocity-V2y= ", V2y*uv0, "km/s"
     print*,"   Velocity-V2z= ", V2z*uv0, "km/s"  
!     print*,""
     print*,"             V2= ", V2*uv0, "km/s"
     print*,""

! ------ Start Universal Variables refining ----------------
!      aprec=a ! previous value of iteration
      sigma2 = (x2 * V2x + y2 * V2y + z2 * V2z)/SQRT(MUc)           ! Dot product of vector RR2,V2
!     
      a = 1D0 /(2D0/r2 - V2**2D0/MUc)          ! a= semiax     
      alfa = -MUc / a   
      print*,"              a= ",a*RT," km" 
      print*," (delta-a)=", ABS(a-aprec)*RT, " km"
      print*,""

      CALL uv (r2,tau3,sigma2,alfa,ff,gg )
      ff3 = ff
      gg3 = gg
      print*,'   (f3,g3)= ', ff3,gg3
!      print*,""

      CALL uv (r2,tau1,sigma2,alfa,ff,gg )
      ff1 = ff
      gg1 = -gg
      print*,'   (f1,g1)= ', ff1,gg1
!      print*,""
    
      c1 = gg3 / (ff1 * gg3 - ff3 * gg1)
      c3 =-gg1 / (ff1 * gg3 - ff3 * gg1)
      print*,'   (c1,c3)= ', c1,c3
      print*,""

!     -------- check cycle -------
        IF (ABS(a-aprec)*RT > 1d-5) then   ! refination down to cm's !!
          goto 1234
        else
          goto 9999
        end if

! -------------------
9999  continue
     
      CALL pveq (x2*RT,y2*RT,z2*RT,V2x*uv0,V2y*uv0,V2z*uv0)  ! write directly keplerian elements of orbit 
	  WRITE(*,*)"*********************************************************"
      WRITE(*,*)"             Press any Key to EXIT                       "
      WRITE(*,*)"*********************************************************"

     READ *, SIGN

 end program orbitsei



