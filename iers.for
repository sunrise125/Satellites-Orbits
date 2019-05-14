*-----------------------------------------------------------------------
      subroutine iers_calc(UTC,DATE,FINAL,XP,YP,DUT1,LOD,dX,dY)

C
C This subroutine detects POLAR MOTION (XP, YP, UT1-UTC, dX and DY from
C IERS database "iers2000A all" (1973-2012) and passes the values
C for the subsequent interpolation and computation of diurnal gravitational  
C and tidal ocanic effects.                
C                    coded by Aldo Nicola (revised Set.2012)  

C choose between A or B bulletin IERS final_all (IAU2000) Standard EOP Data
C Bulletin A , between lines 38 - 46


      implicit none

      DOUBLE PRECISION  XP,YP,DUT1,X1,Y1,RJDINT,UTC,MJD        
      DOUBLE PRECISION YEAR,X(7),Y(7),UT1(7),rjd_int,RJD(7)
      DOUBLE PRECISION JD,DATE,Xx,Yy,DUT,x_int,y_int,ut1_int
      DOUBLE PRECISION dX,dY,DDX(7),DDY(7),DEX,DEY,lo_d
      DOUBLE PRECISION Ra,Rb,Rc,Rd,Re,Rf,DDT(7),DET,LOD
      INTEGER :: A,I,A1,N,J,FINAL                        
      INTEGER :: aa,bb,cc,ios,loop
      character(len=1):: B,B1,dd,ee,ff
      INTEGER E(7)
      INTEGER F(7)
      INTEGER G(7)

      
      A1 = int(3)         
      A=int(DATE - A1)
      FINAL = FINAL + 3
C     Open database  "iers_finals_2000A" 
      open (unit=13, file="finals-iers.txt", status="old",         
     .  action="read",form="formatted",position="rewind")

C     Open database  "iers1962_now" 
C      open (unit=13, file="iers1962_now.txt", status="old",         
C     .  action="read",form="formatted",position="rewind")


      do loop = 41684,FINAL          


C------ Use iers1962_now (IAU2000) EOP08 C04 Data- 
C      read (unit=13,fmt="(3(i4),i7,2(f11.6),2(f12.7),2(f11.6))",
C     .  iostat=ios)aa,bb,cc,MJD,Xx,Yy,DUT,lod,dX,dY
      
 
      
C------ Use IERS final_all (IAU2000) Standard EOP Data- from Bulletin A
      read (unit=13,fmt="(i2,i2,i2,x,f8.2,x,a1,x,f9.6,f9.6,x,f9.6,f9.6,  
     .  2x,a1,f10.7,a12,f6.4,a14,f6.3,a13,f6.3)",iostat=ios)aa,bb,cc,
     .  MJD,B,Xx,X1,Yy,Y1,B1,DUT,dd,lo_d,ff,dX,ee,dY     
      

C------ Use IERS final_all (IAU2000) Standard EOP Data- from Bulletin B                                     
C      read (unit=13,fmt="(i6,x,f8.2,120x,f9.6,x,f9.6,x,f10.7,3x,f7.3,  
C     .  4x,f6.3)",iostat=ios)aa,MJD,Xx,Yy,DUT,dX,dY

      
      
      if (ios==0) then 
       continue
      else
       exit
      end if
      
      do I = 1,7              !  reads  7  values of each parameter           
       if(loop == A + I) then 
       
        E(I)= aa
        F(I)= bb  
        G(I)= cc
        RJD(I)= MJD
        X(I)= Xx 
        Y(I)= Yy
        UT1(I)= DUT
        DDT(I) = lo_d
        DDX(I)= dX
        DDY(I)= dY  
       
       end if
      end do
      end do
      
      rewind (unit=13)
      
C    " N ,NUMBER OF PAIRS FOR INTERPOLATION ( X, Y)"
      N = 7

C    "INTERPOLATION FOR UTC (MJD)
      rjd_int = UTC
      RJD(N) = RJD(I)
      X(N) = X(I)
      Y(N) = Y(I)
      UT1(N) = UT1(I)
      DDT(N) = DDT(I)
      DDX(N) = DDX(I)
      DDY(N) = DDY(I)
C      DDT = 0D0
C     
C    Calling Sub INTERP,LAGINT,PMUT1_OCEANS,PM_GRAVI
C
      
      CALL INTERP(RJD,X,Y,UT1,DDT,DDX,DDY,N,rjd_int,x_int,y_int,ut1_int, 
     .  DET,DEX,DEY) 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C                                                                           !
C  Recoding values of Polar motion (XP , YP), (UT1-UTC), LOD and (dX , dY)  !                                                            !
C                                                                           ! 
C!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      XP = x_int
      YP = y_int
      DUT1 = ut1_int
      dX = DEX
      dY = DEY
      LOD = DET 
!      close (13)            ! added Mar.30,2013  (Giuseppe)
      end 
C-----------------------------------------------------------------------------------------------------------

C----------------------------------------------------------------
      SUBROUTINE INTERP(RJD,X,Y,UT1,DDT,DDX,DDY,N,rjd_int,  
     .  x_int,y_int,ut1_int,DET,DEX,DEY)
C
C     This subroutine takes a series of x, y, and UT1-UTC values
C     and interpolates them to an epoch of choice. This routine
C     assumes that the values of x and y are in seconds of
C    arc and that UT1-UTC is in seconds of time. At least
C     one point before and one point after the epoch of the
C     interpolation point are necessary in order for the
C     interpolation scheme to work. 
C
C     parameters are :
C     RJD     - array of the epochs of data (given in mjd)
C     X       - array of x polar motion (arcsec)
C     Y       - array of y polar motion (arcsec)
C     UT1     - array of UT1-UTC (sec)
C     n       - number of points in arrays
C     rjd_int - epoch for the interpolated value
C     x_int   - interpolated value of x
C     y_int   - interpolated value of y
C     ut1_int - interpolated value of ut1-utc
C
C     CALLED SUBROUTINE : LAGINT (Lagrange interpolation)
C                         PMUT1_OCEANS (Diurnal and semidiurnal oceanic effects)
C                         PM_GRAVI (Diurnal and semidiurnal lunisolar effects)
C
C      coded by Ch. BIZOUARD (Observatoire de Paris) : November 2002
C                                          Corrected : September 2007   
  
      implicit none
      integer I,N
      DOUBLE PRECISION RJD(N), X(N), Y(N), UT1(N),DDX(N),DDY(N),DDT(N),
     . rjd_int, x_int, y_int, ut1_int, DEX, DEY,DET,lod,
     . cor_x, cor_y, cor_ut1, cor_lod
      
       
      CALL LAGINT (RJD,X,n,rjd_int,x_int)
  
      CALL LAGINT (RJD,Y,n,rjd_int,y_int)
      
      CALL LAGINT (RJD,UT1,n,rjd_int,ut1_int)

      CALL LAGINT (RJD,DDX,n,rjd_int,DEX)

      CALL LAGINT (RJD,DDY,n,rjd_int,DEY)

      CALL LAGINT (RJD,DDT,n,rjd_int,DET)
      
C --------------
C Oceanic effect      
C --------------
   
      CALL PMUT1_OCEANS (rjd_int,cor_x,cor_y,cor_ut1,cor_lod)

      x_int = x_int + cor_x
      y_int = y_int + cor_y
      ut1_int = ut1_int + cor_ut1
      lod = lod + cor_lod

C Lunisolar effect 
      CALL PM_GRAVI (rjd_int,cor_x,cor_y)
      
      x_int   = x_int + cor_x
      y_int   = y_int + cor_y
     
      RETURN

      END
C
C----------------------------------------------------------------
C



      SUBROUTINE LAGINT (X,Y,n,xint,yout)
C 
C     This subroutine performs lagrangian interpolation
C     within a set of (X,Y) pairs to give the y
C     value corresponding to xint. This program uses a
C     window of 4 data points to perform the interpolation.
C     if the window size needs to be changed, this can be
C     done by changing the indices in the do loops for
C     variables m and j.
C
C     PARAMETERS ARE :
C     X     - array of values of the independent variable
C     Y     - array of function values corresponding to x
C     n     - number of points
C     xint  - the x-value for which estimate of y is desired
C     yout  - the y value returned to caller

      implicit none
 
      DOUBLE PRECISION X(n),Y(n),xint,yout,term
      INTEGER i,j,k,m,n

      yout = 0.d0
      do  i = 1,n-1
        if ( xint .ge. X(i) .and. xint .lt. X(i+1) ) k = i
      enddo
    
      if ( k .lt. 2 ) k = 2
      if ( k .gt. n-2 ) k = n-2
   
      do m = k-1,k+2
        term = y(m)
        do  j = k-1,k+2
          if ( m .ne. j ) then
            term = term * (xint - X(j))/(X(m) - X(j))
          end if
        enddo 
        yout = yout + term
      enddo
   
      return
      end

C----------------------------------------------------------------
      SUBROUTINE PMUT1_OCEANS (rjd,cor_x,cor_y,cor_ut1,cor_lod)
C
C    This subroutine provides, in time domain, the diurnal/subdiurnal
C    tidal effets on polar motion ("), UT1 (s) and LOD (s). The tidal terms,
C    listed in the program above, have been extracted from the procedure   
C    ortho_eop.f coed by Eanes in 1997.
C    
C    N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
C
C    These corrections should be added to "average"
C    EOP values to get estimates of the instantaneous values.
C
C     PARAMETERS ARE :
C     rjd      - epoch of interest given in mjd
C     cor_x    - tidal correction in x (sec. of arc)
C     cor_y    - tidal correction in y (sec. of arc)
C     cor_ut1  - tidal correction in UT1-UTC (sec. of time)
C     cor_lod  - tidal correction in length of day (sec. of time)
C
C     coded by Ch. Bizouard (2002), initially coded by McCarthy and 
C     D.Gambis(1997) for the 8 prominent tidal waves.  
      
      IMPLICIT NONE
      
      INTEGER nlines
      PARAMETER(nlines=71)
      DOUBLE PRECISION ARG(6),    ! Array of the tidal arguments   
     .                 DARG(6)    ! Array of their time derivative 
      
      REAL*4 XCOS(nlines),XSIN(nlines),
     .YCOS(nlines),YSIN(nlines),UTCOS(nlines),UTSIN(nlines)
     
      DOUBLE PRECISION t,ag,dag,rjd,halfpi,secrad,
     .       cor_x,cor_y,cor_ut1,cor_lod
      INTEGER NARG(nlines,6),i,j
      
      halfpi = 1.5707963267948966d0
      secrad=2.d0*halfpi/(180.d0*3600.d0)	

c  Oceanic tidal terms present in x (microas),y(microas),ut1(microseconds)       
c  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
       data( 
     & NARG(j,1),NARG(j,2),NARG(j,3),NARG(j,4),NARG(j,5),NARG(j,6),
     & XSIN(j),XCOS(j),YSIN(j),YCOS(j),UTSIN(j),UTCOS(j),j=1,nlines)/
     &1,-1, 0,-2,-2,-2,  -0.05,   0.94,  -0.94,  -0.05,  0.396, -0.078,
     &1,-2, 0,-2, 0,-1,   0.06,   0.64,  -0.64,   0.06,  0.195, -0.059,
     &1,-2, 0,-2, 0,-2,   0.30,   3.42,  -3.42,   0.30,  1.034, -0.314,
     &1, 0, 0,-2,-2,-1,   0.08,   0.78,  -0.78,   0.08,  0.224, -0.073,
     &1, 0, 0,-2,-2,-2,   0.46,   4.15,  -4.15,   0.45,  1.187, -0.387,
     &1,-1, 0,-2, 0,-1,   1.19,   4.96,  -4.96,   1.19,  0.966, -0.474,
     &1,-1, 0,-2, 0,-2,   6.24,  26.31, -26.31,   6.23,  5.118, -2.499,
     &1, 1, 0,-2,-2,-1,   0.24,   0.94,  -0.94,   0.24,  0.172, -0.090,
     &1, 1, 0,-2,-2,-2,   1.28,   4.99,  -4.99,   1.28,  0.911, -0.475,
     &1, 0, 0,-2, 0, 0,  -0.28,  -0.77,   0.77,  -0.28, -0.093,  0.070,
     &1, 0, 0,-2, 0,-1,   9.22,  25.06, -25.06,   9.22,  3.025, -2.280,
     &1, 0, 0,-2, 0,-2,  48.82, 132.91,-132.90,  48.82, 16.020,-12.069,
     &1,-2, 0, 0, 0, 0,  -0.32,  -0.86,   0.86,  -0.32, -0.103,  0.078,
     &1, 0, 0, 0,-2, 0,  -0.66,  -1.72,   1.72,  -0.66, -0.194,  0.154,
     &1,-1, 0,-2, 2,-2,  -0.42,  -0.92,   0.92,  -0.42, -0.083,  0.074,
     &1, 1, 0,-2, 0,-1,  -0.30,  -0.64,   0.64,  -0.30, -0.057,  0.050,
     &1, 1, 0,-2, 0,-2,  -1.61,  -3.46,   3.46,  -1.61, -0.308,  0.271,
     &1,-1, 0, 0, 0, 0,  -4.48,  -9.61,   9.61,  -4.48, -0.856,  0.751,
     &1,-1, 0, 0, 0,-1,  -0.90,  -1.93,   1.93,  -0.90, -0.172,  0.151,
     &1, 1, 0, 0,-2, 0,  -0.86,  -1.81,   1.81,  -0.86, -0.161,  0.137,
     &1, 0,-1,-2, 2,-2,   1.54,   3.03,  -3.03,   1.54,  0.315, -0.189,
     &1, 0, 0,-2, 2,-1,  -0.29,  -0.58,   0.58,  -0.29, -0.062,  0.035,
     &1, 0, 0,-2, 2,-2,  26.13,  51.25, -51.25,  26.13,  5.512, -3.095,
     &1, 0, 1,-2, 2,-2,  -0.22,  -0.42,   0.42,  -0.22, -0.047,  0.025,
     &1, 0,-1, 0, 0, 0,  -0.61,  -1.20,   1.20,  -0.61, -0.134,  0.070,
     &1, 0, 0, 0, 0, 1,   1.54,   3.00,  -3.00,   1.54,  0.348, -0.171,
     &1, 0, 0, 0, 0, 0, -77.48,-151.74, 151.74, -77.48,-17.620,  8.548,
     &1, 0, 0, 0, 0,-1, -10.52, -20.56,  20.56, -10.52, -2.392,  1.159,
     &1, 0, 0, 0, 0,-2,   0.23,   0.44,  -0.44,   0.23,  0.052, -0.025,
     &1, 0, 1, 0, 0, 0,  -0.61,  -1.19,   1.19,  -0.61, -0.144,  0.065,
     &1, 0, 0, 2,-2, 2,  -1.09,  -2.11,   2.11,  -1.09, -0.267,  0.111,
     &1,-1, 0, 0, 2, 0,  -0.69,  -1.43,   1.43,  -0.69, -0.288,  0.043,
     &1, 1, 0, 0, 0, 0,  -3.46,  -7.28,   7.28,  -3.46, -1.610,  0.187,
     &1, 1, 0, 0, 0,-1,  -0.69,  -1.44,   1.44,  -0.69, -0.320,  0.037,
     &1, 0, 0, 0, 2, 0,  -0.37,  -1.06,   1.06,  -0.37, -0.407, -0.005,
     &1, 2, 0, 0, 0, 0,  -0.17,  -0.51,   0.51,  -0.17, -0.213, -0.005,
     &1, 0, 0, 2, 0, 2,  -1.10,  -3.42,   3.42,  -1.09, -1.436, -0.037,
     &1, 0, 0, 2, 0, 1,  -0.70,  -2.19,   2.19,  -0.70, -0.921, -0.023,
     &1, 0, 0, 2, 0, 0,  -0.15,  -0.46,   0.46,  -0.15, -0.193, -0.005,
     &1, 1, 0, 2, 0, 2,  -0.03,  -0.59,   0.59,  -0.03, -0.396, -0.024,
     &1, 1, 0, 2, 0, 1,  -0.02,  -0.38,   0.38,  -0.02, -0.253, -0.015,
     &2,-3, 0,-2, 0,-2,  -0.49,  -0.04,   0.63,   0.24, -0.089, -0.011,
     &2,-1, 0,-2,-2,-2,  -1.33,  -0.17,   1.53,   0.68, -0.224, -0.032,
     &2,-2, 0,-2, 0,-2,  -6.08,  -1.61,   3.13,   3.35, -0.637, -0.177,
     &2, 0, 0,-2,-2,-2,  -7.59,  -2.05,   3.44,   4.23, -0.745, -0.222,
     &2, 0, 1,-2,-2,-2,  -0.52,  -0.14,   0.22,   0.29, -0.049, -0.015,
     &2,-1,-1,-2, 0,-2,   0.47,   0.11,  -0.10,  -0.27,  0.033,  0.013,
     &2,-1, 0,-2, 0,-1,   2.12,   0.49,  -0.41,  -1.23,  0.141,  0.058,
     &2,-1, 0,-2, 0,-2, -56.87, -12.93,  11.15,  32.88, -3.795, -1.556,
     &2,-1, 1,-2, 0,-2,  -0.54,  -0.12,   0.10,   0.31, -0.035, -0.015,
     &2, 1, 0,-2,-2,-2, -11.01,  -2.40,   1.89,   6.41, -0.698, -0.298,
     &2, 1, 1,-2,-2,-2,  -0.51,  -0.11,   0.08,   0.30, -0.032, -0.014,
     &2,-2, 0,-2, 2,-2,   0.98,   0.11,  -0.11,  -0.58,  0.050,  0.022,
     &2, 0,-1,-2, 0,-2,   1.13,   0.11,  -0.13,  -0.67,  0.056,  0.025,
     &2, 0, 0,-2, 0,-1,  12.32,   1.00,  -1.41,  -7.31,  0.605,  0.266,
     &2, 0, 0,-2, 0,-2,-330.15, -26.96,  37.58, 195.92,-16.195, -7.140,
     &2, 0, 1,-2, 0,-2,  -1.01,  -0.07,   0.11,   0.60, -0.049, -0.021,
     &2,-1, 0,-2, 2,-2,   2.47,  -0.28,  -0.44,  -1.48,  0.111,  0.034,
     &2, 1, 0,-2, 0,-2,   9.40,  -1.44,  -1.88,  -5.65,  0.425,  0.117,
     &2,-1, 0, 0, 0, 0,  -2.35,   0.37,   0.47,   1.41, -0.106, -0.029,
     &2,-1, 0, 0, 0,-1,  -1.04,   0.17,   0.21,   0.62, -0.047, -0.013,
     &2, 0,-1,-2, 2,-2,  -8.51,   3.50,   3.29,   5.11, -0.437, -0.019,
     &2, 0, 0,-2, 2,-2,-144.13,  63.56,  59.23,  86.56, -7.547, -0.159,
     &2, 0, 1,-2, 2,-2,   1.19,  -0.56,  -0.52,  -0.72,  0.064,  0.000,
     &2, 0, 0, 0, 0, 1,   0.49,  -0.25,  -0.23,  -0.29,  0.027, -0.001,
     &2, 0, 0, 0, 0, 0, -38.48,  19.14,  17.72,  23.11, -2.104,  0.041,
     &2, 0, 0, 0, 0,-1, -11.44,   5.75,   5.32,   6.87, -0.627,  0.015,
     &2, 0, 0, 0, 0,-2,  -1.24,   0.63,   0.58,   0.75, -0.068,  0.002,
     &2, 1, 0, 0, 0, 0,  -1.77,   1.79,   1.71,   1.04, -0.146,  0.037,
     &2, 1, 0, 0, 0,-1,  -0.77,   0.78,   0.75,   0.45, -0.064,  0.017,
     &2, 0, 0, 2, 0, 2,  -0.33,   0.62,   0.65,   0.19, -0.049,  0.018/

      T = (rjd - 51544.5D0)/36525.0D0  ! julian century

C Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
C et leur derivee temporelle 

      ARG(1) = (67310.54841d0 +
     .        (876600d0*3600d0 + 8640184.812866d0)*T +
     .         0.093104d0*T**2 -
     .         6.2d-6*T**3)*15.0d0 + 648000.0d0
      ARG(1)=dmod(ARG(1),1296000d0)*secrad 
   
      DARG(1) = (876600d0*3600d0 + 8640184.812866d0 
     .         + 2.d0 * 0.093104d0 * T - 3.d0 * 6.2d-6*T**2)*15.d0
      DARG(1) = DARG(1)* secrad / 36525.0D0   ! rad/day


      ARG(2) = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2
     .  + 1717915923.2178d0*T + 485868.249036d0
      ARG(2) = DMOD(ARG(2),1296000d0)*secrad
      
      DARG(2) = -4.d0*0.00024470d0*T**3 + 3.d0*0.051635d0*T**2 
     .  + 2.d0*31.8792d0*T + 1717915923.2178d0 
      DARG(2) = DARG(2)* secrad / 36525.0D0   ! rad/day

      ARG(3) = -0.00001149d0*T**4 - 0.000136d0*T**3
     .  -  0.5532d0*T**2
     .  + 129596581.0481d0*T + 1287104.79305d0
      ARG(3) = DMOD(ARG(3),1296000d0)*secrad

      DARG(3) = -4.D0*0.00001149d0*T**3 - 3.d0*0.000136d0*T**2
     .  -  2.D0*0.5532d0*T + 129596581.0481d0
      DARG(3) = DARG(3)* secrad / 36525.0D0   ! rad/day
          
      ARG(4) = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2
     .  + 1739527262.8478d0*T + 335779.526232d0
      ARG(4) = DMOD(ARG(4),1296000d0)*secrad

      DARG(4) = 4.d0*0.00000417d0*T**3 - 3.d0*0.001037d0*T**2 
     .- 2.d0 * 12.7512d0*T + 1739527262.8478d0 
      DARG(4) = DARG(4)* secrad / 36525.0D0   ! rad/day
    
      ARG(5) = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2
     .  + 1602961601.2090d0*T + 1072260.70369d0
      ARG(5) = DMOD(ARG(5),1296000d0)*secrad

      DARG(5) = -4.d0*0.00003169d0*T**3 + 3.d0*0.006593d0*T**2
     . - 2.d0 * 6.3706d0*T + 1602961601.2090d0
      DARG(5) = DARG(5)* secrad / 36525.0D0   ! rad/day

      ARG(6) = -0.00005939d0*T**4 + 0.007702d0*T**3
     .  + 7.4722d0*T**2
     .  - 6962890.2665d0*T + 450160.398036d0
      ARG(6) = DMOD(ARG(6),1296000d0)*secrad

      DARG(6) = -4.d0*0.00005939d0*T**3 + 3.d0 * 0.007702d0*T**2
     .  + 2.d0 * 7.4722d0*T - 6962890.2665d0
      DARG(6) = DARG(6)* secrad / 36525.0D0   ! rad/day

C CORRECTIONS

	cor_x  = 0.d0
	cor_y  = 0.d0
	cor_ut1= 0.d0
	cor_lod= 0.d0

 	do j=1,nlines
 	
	ag  = 0.d0
 	dag = 0.d0
		do i=1,6
 		ag  = ag  + dble(narg(j,i))*ARG(i)
 		dag = dag + dble(narg(j,i))*DARG(i)
		enddo
	ag=dmod(ag,4.d0*halfpi)

        cor_x= cor_x + dble(XCOS(j))*dcos(ag) + dble(XSIN(j)) * dsin(ag)
        cor_y= cor_y + dble(YCOS(j))*dcos(ag) + dble(YSIN(j)) * dsin(ag)
        cor_ut1= cor_ut1+dble(UTCOS(j))*dcos(ag)+dble(UTSIN(j))*dsin(ag)
        cor_lod= cor_lod -(-dble(UTCOS(j)) * dsin(ag) 
     &                    + dble(UTSIN(j)) * dcos(ag) ) * dag   	 

        enddo
  
       cor_x   = cor_x * 1.0d-6   ! arcsecond (")
       cor_y   = cor_y * 1.0d-6   ! arcsecond (")
       cor_ut1 = cor_ut1 * 1.0d-6 ! second (s)
       cor_lod = cor_lod * 1.0d-6 ! second (s)
 
      RETURN
      END
      	
C----------------------------------------------------------------
      SUBROUTINE PM_GRAVI (rjd,cor_x,cor_y)
C
C    This subroutine provides, in time domain, the diurnal
C    lunisolar effet on polar motion (")
C    
C    N.B.:  The fundamental lunisolar arguments are those of Simon et al.  
C
C    These corrections should be added to "average"
C    EOP values to get estimates of the instantaneous values.
C
C     PARAMETERS ARE :
C     rjd      - epoch of interest given in mjd
C     cor_x    - tidal correction in x (sec. of arc)
C     cor_y    - tidal correction in y (sec. of arc)
C
C     coded by Ch. Bizouard (2002)
      
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INTEGER nlines
      PARAMETER(nlines=10)
      DOUBLE PRECISION ARG(6)    ! Array of the tidal arguments   
      REAL*4 XCOS(nlines),XSIN(nlines),YCOS(nlines),YSIN(nlines)
      INTEGER NARG(nlines,6)
      
      halfpi = 1.5707963267948966d0
      secrad=2.d0*halfpi/(180.d0*3600.d0)	

c  Diurnal lunisolar tidal terms present in x (microas),y(microas)      
c  NARG(j,6) : Multipliers of GMST+pi and Delaunay arguments. 
	
       data( 
     & NARG(j,1),NARG(j,2),NARG(j,3),NARG(j,4),NARG(j,5),NARG(j,6),
     & XSIN(j),XCOS(j),YSIN(j),YCOS(j),j=1,nlines)/    
     & 1,-1, 0,-2, 0,-1,    -.44,   .25,   -.25,  -.44,
     & 1,-1, 0,-2, 0,-2,   -2.31,  1.32,  -1.32, -2.31,
     & 1, 1, 0,-2,-2,-2,    -.44,   .25,   -.25,  -.44,
     & 1, 0, 0,-2, 0,-1,   -2.14,  1.23,  -1.23, -2.14,
     & 1, 0, 0,-2, 0,-2,  -11.36,  6.52,  -6.52,-11.36,
     & 1,-1, 0, 0, 0, 0,     .84,  -.48,    .48,   .84,
     & 1, 0, 0,-2, 2,-2,   -4.76,  2.73,  -2.73, -4.76,
     & 1, 0, 0, 0, 0, 0,   14.27, -8.19,   8.19, 14.27,
     & 1, 0, 0, 0, 0,-1,    1.93, -1.11,   1.11,  1.93,
     & 1, 1, 0, 0, 0, 0,     .76,  -.43,    .43,   .76/
 
      T = (rjd - 51544.5D0)/36525.0D0  ! julian century

C Arguments in the following order : chi=GMST+pi,l,lp,F,D,Omega
C et leur derivee temporelle 

      ARG(1) = (67310.54841d0 +
     .        (876600d0*3600d0 + 8640184.812866d0)*T +
     .         0.093104d0*T**2 -
     .         6.2d-6*T**3)*15.0d0 + 648000.0d0
      ARG(1)=dmod(ARG(1),1296000d0)*secrad 
   

      ARG(2) = -0.00024470d0*T**4 + 0.051635d0*T**3 + 31.8792d0*T**2
     .  + 1717915923.2178d0*T + 485868.249036d0
      ARG(2) = DMOD(ARG(2),1296000d0)*secrad
      

      ARG(3) = -0.00001149d0*T**4 - 0.000136d0*T**3
     .  -  0.5532d0*T**2
     .  + 129596581.0481d0*T + 1287104.79305d0
      ARG(3) = DMOD(ARG(3),1296000d0)*secrad

          
      ARG(4) = 0.00000417d0*T**4 - 0.001037d0*T**3 - 12.7512d0*T**2
     .  + 1739527262.8478d0*T + 335779.526232d0
      ARG(4) = DMOD(ARG(4),1296000d0)*secrad

    
      ARG(5) = -0.00003169d0*T**4 + 0.006593d0*T**3 - 6.3706d0*T**2
     .  + 1602961601.2090d0*T + 1072260.70369d0
      ARG(5) = DMOD(ARG(5),1296000d0)*secrad

  
      ARG(6) = -0.00005939d0*T**4 + 0.007702d0*T**3
     .  + 7.4722d0*T**2
     .  - 6962890.2665d0*T + 450160.398036d0
      ARG(6) = DMOD(ARG(6),1296000d0)*secrad


C CORRECTIONS

	cor_x  = 0.d0
	cor_y  = 0.d0

 	do j=1,nlines
 	
	ag  = 0.d0
		do i=1,6
 		ag  = ag  + dble(narg(j,i))*ARG(i)
		enddo
	ag=dmod(ag,4.d0*halfpi)

        cor_x =cor_x+dble(XCOS(j))*dcos(ag)+dble(XSIN(j))*dsin(ag)
        cor_y =cor_y+dble(YCOS(j))*dcos(ag)+dble(YSIN(j))*dsin(ag) 

        enddo
  
      cor_x = cor_x * 1.0d-6   ! arcsecond (")
      cor_y = cor_y * 1.0d-6   ! arcsecond (")
 
      RETURN

      END
!-------------------------------------------------------------------------------------------------------


      SUBROUTINE time_hms ( FD, IH, IM, SEC) 

      IMPLICIT NONE
      DOUBLE PRECISION FD,DH,DM,SEC
      INTEGER IH,IM

      IH = INT( FD * 24D0)        ! transform part of a day ( FD) 
      DH = ( FD *24D0) - IH       ! in hours, minuts and seconds
      IM = INT( DH * 60D0)
      DM = ( DH * 60D0) - IM
      SEC = DM * 60D0
     
      END 


*==============================================================================================
      SUBROUTINE GEOCPV ( ELONG, PHI, HEIGHT, GST, POS, VEL, J )
!     This subroutine compute the velocity vector of a terrestrial observer 
!          with respect to the geocenter.       
      
!       Input:  ELONG  = Longitude of geodetic observer in radians ( + EAST)
!               PHI    = Latitude  of geodetic observer in radians ( + NORD)
!               GST    = Greenwich apparent sideral time in radians
!               HEIGHT = HEIGHT of observer in meter 
!       Output  VEL    = Velocity vector of observer (Geocenter equatorial 
!                          rectangular coordinate in AU/Day)
!               POS    = Position vector of observer (Geocenter equatorial 
!                          rectangular coordinate in AU)
!               J      = status:  0 = OK
!                                -1 = illegal case 
!
      IMPLICIT NONE
      
      DOUBLE PRECISION DR2D
      PARAMETER ( DR2D = 57.29577951308232087679815D0 )  ! Radians to degrees

      DOUBLE PRECISION AU
      PARAMETER ( AU = 149597870.691D0)                    ! Astronomical Unit (Km)(TDB) 
   
!     NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY OF EARTH
!     RADIANS/SECOND, FROM IERS CONVENTIONS (2003)      
      DOUBLE PRECISION OMEGA
      PARAMETER (OMEGA = 7.2921150D-5)                   

!     EQUATORIAL RADIUS OF EARTH IN Km, FROM IERS CONVENTIONS (2003)
      DOUBLE PRECISION A
      PARAMETER ( A = 6378.1366D0 )

!     FLATTENING FACTOR OF EARTH, FROM IERS CONVENTIONS (2003)
      DOUBLE PRECISION F
      PARAMETER ( F = 1.D0 / 298.25642D0 )

      DOUBLE PRECISION ELONG, PHI, HEIGHT, GST, VEL(3),POS(3)
      INTEGER J,I

      DOUBLE PRECISION SP, CP, W, D, AC, AS, R, GP(3),GV(3)


!  Functions of geodetic latitude.
      SP = SIN(PHI)
      CP = COS(PHI)
      W = 1D0-F
      W = W*W
      D = CP*CP + W*SP*SP
      IF ( D .GT. 0D0 ) THEN
         AC = A / SQRT(D)
         AS = W * AC

*     Geocentric position vector.
      R = ( AC + HEIGHT/1000D0 ) * CP     
      POS(1) = R * COS(ELONG+GST)
      POS(2) = R * SIN(ELONG+GST)
      POS(3) = ( AS + HEIGHT/1000D0 ) * SP
      

*     Geocentric velocity vector in KM/SEC
      VEL(1) = -OMEGA * R * sin(ELONG+GST)
      VEL(2) =  OMEGA * R * cos(ELONG+GST)
      VEL(3) =  0.D0

*     Success.
         J = 0
      ELSE

*     Fail.
         J = -1
      END IF

*     CONVERT POSITION AND VELOCITY COMPONENTS TO AU AND AU/DAY
      DO 10 I=1,3
      POS(I) = POS(I) / AU
      VEL(I) = VEL(I) / AU * 86400.D0
   10 CONTINUE 

      END
!----------------------------------------------------------------------------------------------

