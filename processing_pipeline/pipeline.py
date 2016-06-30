import numpy as np
from sklearn import gaussian_process
import astropy

from scipy.signal import argrelmin
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
	
def import_idl(filename = '../../../Data/K16/kic126corr_n.sav', num = 400):
	idlfile = rs(filename);
	ident = np.array(idlfile['cont']).astype('int');
	flux = np.array(idlfile['flux']).astype('float64');
	cadence = np.array(idlfile['time']).astype('float64');

	for a in np.unique(ident): #data quarter
		mean = np.average(flux[idlfile['cont'] == a]);
		data = [cadence[idlfile['cont'] == a], flux[idlfile['cont'] == a]];
		print len(data[0]), len(data[1]);
		k =  movavg_final(data, num);
		
		flux[np.where(idlfile['cont'] == a)[0]] -= k.astype('float64');

		#gp2(np.array([cadence[idlfile['cont'] == a], flux[idlfile['cont'] == a]]), block_size = 2000)[2];
		#flux[idlfile['cont'] == a] -= mean;
		flux[idlfile['cont'] == a] /= mean;

	arm = argrelmin(flux)[0];
	arm = arm[flux[arm] < -0.005]
	for u in arm:
		fluxbase = flux[max(0,u - int(num)):min(len(flux), u + int(num))];
		fluxbase_mean = np.average(fluxbase[fluxbase > 0])
		flux[max(0,u - int(num)):min(len(flux), u + int(num))] -= fluxbase_mean;

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
	#print data[1]
	a = pd.rolling_mean(data[1], movavg_period)
	b = pd.rolling_mean(data[1][::-1], movavg_period)[::-1]
	avg = (a + b)/2

	avg[0:movavg_period] = b[0:movavg_period];
	avg[len(avg) - movavg_period:len(avg)] = a[len(avg) - movavg_period:len(avg)]

	return avg;

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

def import_star(star_id = "kplr011904151", cad = "llc"):
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
			q[1] /= np.average(q[1]);
			q[1] -= np.average(q[1]);
                        base_flux = np.append(base_flux, q[1]);
                        base_cadence = np.append(base_cadence, q[0]);

	return base_cadence, base_flux;

def detrend_star(cadence, flux, num = 20):
        arm = argrelmin(flux)[0];
        arm = arm[flux[arm] < -0.005]
        for u in arm:
                fluxbase = flux[max(0,u - int(num)):min(len(flux), u + int(num))];
                fluxbase_mean = np.average(fluxbase[fluxbase > 0])
                flux[max(0,u - int(num)):min(len(flux), u + int(num))] -= fluxbase_mean;

	return cadence, flux

def mean_confidence_interval(data, confidence=0.95):
	"""shasan, http://stackoverflow.com/questions/15033511/compute-a-confidence-interval-from-sample-data """
	a = 1.0*np.array(data)
	n = len(a)
	m, se = np.mean(a), scipy.stats.sem(a)
	h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
	return m, m-h, m+h	
