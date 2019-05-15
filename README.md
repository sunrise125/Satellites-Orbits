## Satellites-Orbits
What this package performs is better detailed in the attached PDF, as shortly condensed below.

### Precise Orbit Determination of Satellites through 3 Observations (Gauss method) 
Since Jan.1,2009 any computation involving Earth gravitational field should be performed by means of the new precession-nutation model and this paper presents how to determine the orbital elements of a satellite being known times and topocentric polar coordinates from a given location. The position of Intermediate Pole is taken from 
IERS and strongly interpolated to have its parameters precise at each time of the three observations. 

The main program <i>orbit6.f90</i> is written in Fortran90 and uses  <b>SOFA</b> subroutines coming from here (http://www.iausofa.org)

### Compilation
**Step 1**: move into the working folder

**Step 2**: be sure that the compiler, i.e. <b>gfortran</b> is accessible from the console

**Step 3a**: gfortran orbit6.f90 1818.for iers.for twosubs.f90 -o orbit6.exe  (for Windows users)

**Step 3b**: gfortran orbit6.f90 1818.for iers.for twosubs.f90 -o orbit6      (for Unix/Linux users)

### Running orbit6.exe
Executable already included in the package, along with two ancillary files .dll to provide functionality under any Windows system.
At startup is asked to input the datafile, i.e. ISSback.dat, then push enter and results appear on the screen, as shown in <i>ISSback-output.txt</i>. Computed orbital elements are resumed underneath
<PRE>
  ---------- Output Satellite Orbital Elements ---------------------------
      Semi-major axis=    6796.5726182190183       km
    1-   Eccentricity=    3.4276705224868965E-003  dimensionless
    2-   Eccentricity=    3.4276705224868965E-003  dimensionless
          Inclination=    51.669493839385176       degs
          Node (RAAN)=    254.71395472212248       degs
  Argument of Perigee=    107.96761154537377       degs
         Mean anomaly=    303.86709907783660       degs
         True anomaly=    303.54018063104337       degs
  ------------------------------------------------------------------------
</PRE>
