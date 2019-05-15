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
