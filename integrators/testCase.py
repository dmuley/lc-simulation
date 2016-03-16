from polarIntegrator3 import *;
import matplotlib.pyplot as plt
import numpy as np;

light_blocked = [];

for pos in np.linspace(-200,200,178):
	n = 51;

	c = [[100., 0., 0.],[10.,pos,0.1], [10., ];
	y = getTangentsIntersections(c);
	m = generateRadiiThetas(n,y[0], y[1], y[2]);
	f = rd2(m, c, opt = 0);
	oh= groupAndIntegrate(bounds = f, num = n);
	light_blocked.append(oh);
	print pos;
	
r = np.array(light_blocked);
plt.plot(np.linspace(-200,200,178),-r/(np.pi * 100.**2));
plt.plot(np.linspace(-200,200,178), np.linspace(-0.01,-0.01, 178));
plt.xlim(-200,200);
plt.show();