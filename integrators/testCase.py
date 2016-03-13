from polarIntegrator2_1 import *;
import matplotlib.pyplot as plt
import numpy as np;

light_blocked = []

for pos in np.linspace(-200,200,321):
	c = [[100., 0., 0.],[10.,pos,0.00]];
	y = getTangentsIntersections(c);
	m = generateRadiiThetas(21,y[0], y[1], y[2]);
	f = rd2(m, c, opt = 0);
	oh= groupAndIntegrate(bounds = f, num = 20);
	light_blocked.append(oh);
	
plt.plot(np.linspace(-200,200,321),-np.array(light_blocked)/(np.pi * 100.**2),'.');
plt.xlim(-200,200);
plt.show();