import numpy as np
from scipy import stats

def pathGenerator(m_planet = 1., e_planet = 0., a_planet = 1., per_planet = 1., m_moon = 0.1, e_moon = 0., a_moon = 0.5, per_moon = 0.1, time = 100):
	# note: the major axis of the orbit will always be in the x-direction

	'''
	TO INCLUDE: Means of accounting for moon's gravitational effects on planet. (include mass ratio between moon and planet). 
	Plugin to transit modeling system (create location constraints and return location of planets/moons as a function of time).
	Option to return either planet, moon, or center-of-mass path.'''
	t = np.linspace(0, time, 100000);
	dx_planet = a_planet * np.cos(2 * np.pi * t/per_planet);
	dy_planet = a_planet * (1 - e_planet)/(1 + e_planet) * np.sin(2 * np.pi * t/per_planet);
	
	dx_moon = a_moon * np.cos(2 * np.pi * t/per_moon);
	dy_moon = a_moon * (1 - e_moon)/(1 + e_moon) * np.sin(2 * np.pi * t/per_moon);
	
	
	x = -((dx_planet + dx_moon * m_planet/(m_planet + m_moon))- a_planet * np.sqrt(1 - ((1 - e_planet)/(1 + e_planet))**2)/2 - a_moon * np.sqrt(1 - ((1 - e_moon)/(1 + e_moon))**2)/2);
	y = (dy_planet + dy_moon * m_planet/(m_planet + m_moon));
	
	xpl = -((dx_planet - dx_moon * m_moon/(m_planet + m_moon))- a_planet * np.sqrt(1 - ((1 - e_planet)/(1 + e_planet))**2)/2 - a_moon * np.sqrt(1 - ((1 - e_moon)/(1 + e_moon))**2)/2);
	ypl = (dy_planet - dy_moon * m_moon/(m_planet + m_moon));
	
	return (np.array([x, y]), np.array([xpl, ypl]));