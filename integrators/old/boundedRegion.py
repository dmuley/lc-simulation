import numpy as np
import scipy
import matplotlib.pyplot as plt
import time

def circle(r,precision,xoffset,yoffset):
	n = np.linspace(0, 2*np.pi, precision);
	# return (x, y) list
	# as well as weights for integration
	return (np.cos(n) * r + xoffset, np.sin(n) * r + yoffset);
	
def bounds(region, *args):
	radius = args[0];
	precision = args[1];
	xoffset = args[2];
	yoffset = args[3];
	
	rad_input = np.linspace(0,radius,precision+1)
	x_values = np.array([])
	y_values = np.array([])
	chain_rule=np.array([])
	
	for a in range(1,len(rad_input)):
		region_evaluation = region(rad_input[a], precision,xoffset,yoffset)
		x_values = np.append(x_values, region_evaluation[0])
		y_values = np.append(y_values, region_evaluation[1])
		chain_rule = np.append(chain_rule, np.ones(precision) * ((rad_input[a]**2 - rad_input[a - 1]**2) * np.pi)/precision)
		
	chain_rule = chain_rule
		
	return (x_values, y_values, chain_rule)
	
def polar_riemann_integral(x, y, chain_rule, radius_large, limb_darkening=0.5):	
	polar_integral = 0.;

	#This is the key part of the integrator, you can change the limb darkening law here to whatever you'd like.	
	addition =  np.nan_to_num(np.sqrt(radius_large**2 - x**2 - y**2) * limb_darkening + radius_large * (1 - limb_darkening));
			
	addition2 = addition[np.where(addition > 0)] * chain_rule[np.where(addition > 0)]
	polar_integral = np.sum(addition2)/2
	
	return polar_integral;
	
def circular_restriction(x, y, chain_rule,xoffset,yoffset,radius):
	pl = np.nan_to_num(np.sqrt((x-xoffset)**2 + (y-yoffset)**2))
	i = np.intersect1d(np.where(pl <= radius)[0], np.where(pl > 0)[0]);
	
	return (x[i], y[i],chain_rule[i])


'''	THIS CODE SHOULD ONLY RUN BY HAND. NOT ORDINARILY.
q0 = bounds(circle, 10., 500, 0, 0);

for ldc in np.linspace(0, 1, 21):
	ya = np.array([])


	ti = time.time()
	star = polar_riemann_integral(q0[0], q0[1], q0[2], 10., ldc);
	for re in np.linspace(-11,0,550):
		cir = bounds(circle, 1., 50., re, 0.)
		t = polar_riemann_integral(cir[0], cir[1], cir[2], 10., ldc)
		ya = np.append(ya, t)
	
	tf = time.time()
	print "Final time: " + str(tf - ti) + " s."
	
	plt.plot(np.linspace(-11, 11, 1100), -np.append(ya, ya[::-1])/star)
	
plt.show()
'''
