# Starlink: huge discrepancy on 64 positions from issued TLE 

Reference is made to file starlink_29may.txt, attached in this folder and downloaded on May.29,2019 from NORAD repository http://celestrak.com/NORAD/elements/starlink.txt

Calculations have been performed by means of Pascal code named SATENG3.PAS (here attached) through SGP Units coming from the public domain software http://celestrak.com/software/vallado-sw.php

The processed position of each sat is shown in 111.txt, and here analyzed with the help of a couple of them. There is a strange alternance <i>Visible/NOT-Visible</i> which is <b>discrepant</b> with the observed objects, i.e. Object A having h= -19.4 degs and Object E h= 48.9, visible.
<PRE>
 Observer Coords. (lat,lon,H)= 37.0328  15.0650  370.0
OBJECT A                
1 44235U 19029A   19147.55729408  .00003298  00000-0  89486-4 0  9992
2 44235  52.9978 158.9310 0001303 181.8828 296.9785 15.40623568  1490
 JD0= 2458631.05729408    dt= 2.54131703   ideep= 0
 2019/5/30.09861 Azim.[deg]= 137.0164   Height[deg]= -19.4264   NOT-Visible
 Topoc.Dist.(`range`) [km]=  5328.071
 (X,Y,Z)={    5786   -3610    -181} -> r_geoc=   6823 km
------------------
OBJECT E                
1 44239U 19029E   19147.55806232  .00002267  00000-0  70150-4 0  9990
2 44239  52.9963 158.9546 0002970 202.5743 258.8860 15.37972393  1497
 JD0= 2458631.05806232    dt= 2.54054879   ideep= 0
 2019/5/30.09861 Azim.[deg]= 198.0158   Height[deg]=  48.9160   Visible
 Topoc.Dist.(`range`) [km]=   587.309
 (X,Y,Z)={    2547   -5070    3794} -> r_geoc=   6825 km
</PRE>

The problem is: <b>why</b> does it happen?
