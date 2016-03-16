# lc-simulation

This is a program that simulates lightcurves for systems with multiple planets and moons, based on the results of a nested-Keplerian simulation. The nested-Keplerian computes the positions of planets and moons with a variety of positions, inclinations, orbital phases, periods, etc., which is then converted into an (x,  y) sky-projected position. This in turn is passed to the light-curve integrator, which computes the light curve.

The processing\_pipeline folder contains an earlier pipeline that works only for a single-planet, single-star case. This code can be run out-of-the-box but is slow and of limited capability. As of now the nested-Keplerian does not directly feed into the light-curve integrator.

In order to use the integration to generate a light curve, please make client code along the lines of testCase.py. I do not recommend using a step size below 50 due to the right-hand-side jittering issue. Nevertheless, both sides are fairly accurate, the LHS to within 0.01% and the RHS to within about 0.1%.

_Please do not look into folders marked "old", "misc", or similar. Much of that code is experimental or brief tried-and-failed snippets which may not even run, or if it does may be too slow, untested or unreliable to be useful._

Todo list:

+ Seamless transition from nested-Keplerian to integrator, optimize selections of planets and stars to ensure minimal computation time. Treat the multiple-star case.
+ Fixing excess jittering on the right-hand-side (-pi/2 < theta < pi/2) of a star.
+ Test integrator with various limb-darkening laws (currently only constant one is used for testing purposes).
+ Make further optimizations and accuracy adjustments in the code (use numpy wherever possible, unroll excessively deep/long loops, etc.)
+ Comment code in-depth to make usage and understanding easier.
