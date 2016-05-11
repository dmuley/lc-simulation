import numpy as np
from sklearn import gaussian_process
import astropy

from scipy.io import readsav as rs
from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
import time
import os

def import_data(filename):
	hdulist = fits.open(filename);
	sap_flux = hdulist[1].data["PDCSAP_FLUX"];
	cadence = hdulist[1].data["CADENCENO"];
	
	sf = sap_flux[np.where(np.logical_not(np.isnan(sap_flux)))];
	ca = cadence[np.where(np.logical_not(np.isnan(sap_flux)))];
	
	return np.array([ca, sf]);
	
def import_idl(filename = '../../../Data/K16/kic126corr_n.sav'):
	idlfile = rs(filename);
	flux = idlfile['flux'];
	cadence = idlfile['time'];

	for a in np.unique(idlfile['cont']): #data quarter
		flux[idlfile['cont'] == a] /= np.average(flux[idlfile['cont'] == a]);

	return cadence, flux;

def data_folder(data, period, cadence_period=58.84876 * 30):
	modulus = (period * 3600. * 24./cadence_period);
	print modulus;
	new = data[0] % modulus;
	l = np.argsort(new);
	q = new[l];
	r = data[1][l];
	
	return np.array([q, r]);
	
def gp(data, nug=0.01):
	X = np.atleast_2d(data[0]).T
	y = data[1].ravel()
	x = np.atleast_2d(data[0]).T
	gp = gaussian_process.GaussianProcess(theta0=1e-2, thetaL=1e-4, thetaU=0.5e-1, nugget = nug)
	gp.fit(X, y)  
	y_pred, sigma2_pred = gp.predict(x, eval_MSE=True)
	
	return np.array([data[0], y_pred]);
	
def gp2(data, block_size = 100, nugget = 0.005):

	c = data[0]
	s = data[1]
	s_2 = np.array_split(s, len(s)/block_size + 1)
	c_2 = np.array_split(c, len(s)/block_size + 1)
	
	sapflux_pred = []
	
	nug = nugget;
	for a in range(0,len(s_2)):
	
		t0 = time.time()
		X = np.atleast_2d(c_2[a]).T
		y = np.atleast_2d(s_2[a]).T
	
		gproc = gaussian_process.GaussianProcess(theta0=0.01, thetaL=1e-4, thetaU=1e-1,nugget=nug)
	
		
		gproc.fit(X, y)
		y_pred, sigma2_pred = gproc.predict(X, eval_MSE=True)
		sapflux_pred.extend(y_pred.ravel())
		t1 = time.time()
		print t1-t0
	
	return np.array([c, s, np.array(sapflux_pred)])

	
def movavg_regression(data, movavg_period):
	a = pd.rolling_mean(data[1], movavg_period)
	b = pd.rolling_mean(data[1][::-1], movavg_period)[::-1]
	avg = (a + b)/2
	
	return np.array([data[0][movavg_period:len(avg) - movavg_period], data[1][movavg_period:len(avg) - movavg_period] - avg[movavg_period:len(avg) - movavg_period]]);
	
def ema(data, hl=2):
	a = pd.ewma(data[1], halflife=hl)
	b = pd.ewma(data[1][::-1], halflife=hl)[::-1]
	avg = (a + b)/2
	
	return np.array([data[0], (data[1] - avg)]);
	
def movavg_final(data, movavg_period):
	a = pd.rolling_mean(data[1], movavg_period)
	b = pd.rolling_mean(data[1][::-1], movavg_period)[::-1]
	avg = (a + b)/2

	return np.array([data[0][movavg_period:len(avg) - movavg_period], avg[movavg_period:len(avg) - movavg_period]])

def full_pipeline(star_id = "kplr011904151", cad = "llc", per = 0.837495):
	#SAMPLE CODE THAT RUNS ALL OF THIS.
	#Do not run ordinarily, only by hand. This is a specific testing case.
	base_flux = np.array([])
	base_cadence = np.array([])
	if cad == "llc":
		cad_factor = 30.
	else:
		cad_factor = 1.
	for a in os.listdir("."):
		if star_id and cad in a:
			print a;
		
			q = import_data(a)
			r = movavg_regression(q, int(np.sqrt(len(q[0]))));
			#r = gp2(q, 500)
			#r = q
		
			base_flux = np.append(base_flux, r[1]);
			base_cadence = np.append(base_cadence, r[0]);
		
	s = np.array([base_cadence, base_flux]);
	t = data_folder(s, per, 58.84876 * cad_factor);
	u = movavg_final(t, 2);
	#plt.plot(t[0], t[1], '.');
	#plt.plot(u[0], u[1], '.');
	
	#plt.show();
	
	return t

		
	
