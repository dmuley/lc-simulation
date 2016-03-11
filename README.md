# lc-simulation

This is a program that simulates lightcurves for systems with multiple planets and moons, based on the results of a nested-Keplerian simulation. The nested-Keplerian computes the positions of planets and moons with a variety of positions, inclinations, orbital phases, periods, etc., which is then converted into an (x,  y) sky-projected position. This in turn is passed to the light-curve integrator, which computes the light curve.

The processing\_pipeline folder contains an earlier pipeline that works only for a single-planet, single-star case. This code can be run out-of-the-box but is slow and of limited capability. As of now the nested-Keplerian does not directly feed into the light-curve integrator.

Todo list:

+ Seamless transition from nested-Keplerian to integrator, optimize selections of planets and stars to ensure minimal computation time. Treat the multiple-star case.
+ Test integrator with various limb-darkening laws (currently only constant one is used for testing purposes).
+ Make further optimizations and accuracy adjustments in the code (use numpy wherever possible, unroll excessively deep/long loops, etc.)
+ Fully and in-depth comment code.
