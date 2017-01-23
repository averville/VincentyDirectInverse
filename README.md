# VincentyDirectInverse
Vincenty Direct and Inverse geodetic computations
-------------------------------------------------
Thaddeus Vincenty has put a final point and solved this complex problem in 1975
https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf

The direct problem:
-------------------
From a geodetic position known in latitude and longitude
Given a geodetic azimuth (bearing) and distance to another point
Compute the geodetic latitude and longitude of the second point

The inverse problem:
--------------------
From two geodetic positions known in latitude and longitude
Compute the direct and inverse azimuths (bearings) and distance

Notes:
------
Seems easy, but the computation has to be made on the surface of an
ellipsoid, and be accurate to the positional fraction of a millimeter,
even at very long distances. The direct computation requires a small
iteration to converge with sufficient accuracy.

I have made my "end of studies" thesis on Vincenty's method, comparing it
at the same time to another method from Hradilek.

Originally written for the TI-58 calculator in 1978, I needed them on
a flight navigation project, so I decided to code them in Python.

Here we go. Simple to use. Ellipsoid parameters can be changed easily
within the code or transferred if needed through the provided parameters.
The formulas have been tested and work fine with a few geodetic positions,
hoping they also work in all World areas, but did not test all cases.
