import numpy as np

def var_analysis(x ,y , period_len):
	n = 0;
	final_array = np.array([])
	for i in range(0,len(x)):
		if x[i] < period_len:
			final_array = np.append(final_array, np.nan);
			n = i;
		else:
			processing_array_indices = np.where((x[n:i] - x[i] > -period_len));
			pa = x[processing_array_indices];
			q = np.var(pa);
			final_array = np.append(final_array, q);
			print (q, i);
	
	return final_array;
	
'''

Working on a variance-analysis model of detecting extrasolar planets. Right now we use period folding to detect planets, which I am working on as well. However this doesn't work very well for longer-period planets which are likely to host life, because there are too few periods for the folding's averaging effect to work well. Also this works poorly for nearly all stars except the least variable ones.

An alternate technique may be to analyze the variance of the stellar data (perhaps using an exponential moving average, or similar model) and find points where the variance has been "pinched". These indicate that some part of the variation has been obscured (i. e. the star has been blocked in part by a planet.)

''