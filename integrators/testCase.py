from polarIntegrator3 import *;
import matplotlib.pyplot as plt
import numpy as np;

light_blocked = [];

for pos in np.linspace(-200,200,178):
	n = 71;

	c = [[100., 0., 0.],[10.,pos,0.1]];
	y = getTangentsIntersections(c);
	m = generateRadiiThetas(n,y[0], y[1], y[2]);
	f = rd2(m, c, opt = 0);
	oh= groupAndIntegrate(bounds = f, num = n, star_rad = 100., ld_coeff = 0.);
	light_blocked.append(oh);
	print pos;
	
r = np.array(light_blocked);
plt.plot(np.linspace(-200,200,178),-r, '.');
plt.plot(np.linspace(-200,200,178), np.linspace(-0.01,-0.01, 178));
plt.xlim(-200,200);
plt.show();