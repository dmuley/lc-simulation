# Kepler Generalized Light Curve Simulation

The Kepler Generalized Light Curve Simulation computes light output with respect to time for systems with multiple planets, stars and moons, based on the results of a nested-Keplerian simulation. The nested-Keplerian computes the positions of planets and moons with a variety of positions, inclinations, orbital phases, periods, etc., which is then converted into an (x,  y) sky-projected position. This in turn is passed to the light-curve integrator, which calculates the amount of starlight blocked at each timestep and so constructs the lightcurve.

The processing\_pipeline folder contains an earlier pipeline that works only for a single-planet, single-star case. This code can be run out-of-the-box but is slow and of limited capability. As of now, the nested-Keplerian code has been updated to interface directly with the light-curve integrator, meaning that it is possible to create usable and physically accurate light curves. In order to do this, please write a script along the lines of kepler16.py (the calls to general\_transit.py can remain the same). *I do not recommend using an integration step number below 30* due to jittering on the right-hand side, although lower numbers can be used for sanity checks.

##Features

+ Seamless transition from nested-Keplerian to integrator. Optimize selections of planets and stars for integration (does not calculate light curves of planets) to ensure minimal computation time. 
+ Low jittering on the right-hand-side (-pi/2 < theta < pi/2) of a star.
+ Works with a wide range of limb-darkening laws.
+ Retrograde orbits (add pi radians to inclination and argument of periastron).

##Todo

+ Make further optimizations and accuracy adjustments in the code (use numpy wherever possible, unroll excessively deep/long loops, etc.)
+ Comment code in-depth to make usage and understanding easier.
+ Process jagged OrbitingSystem arrays without using dummy planets or moons.
+ Produce professional-quality visualizations to demonstrate working of program.
+ Test against detrended light-curves from _Kepler_ space telescope in order to ensure veracity of algorithm.

For a more complete and up-to-date list of issues, please look at the Issues tab, and add in any issues that you notice.
