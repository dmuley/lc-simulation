import numpy as np
import matplotlib.pyplot as plt

#Arguments of this function will be each [radius, X, Y];
def circleIntersectionCalculator(list_of_circles):
	left_coords = [];
	right_coords = [];
	y_coords = [];
	for a in list_of_circles:
		#get leftmost and rightmost points for the limits of integration
		left_coords = np.append(left_coords, a[1] - a[0]);
		right_coords= np.append(right_coords, a[1] + a[0]);
		y_coords = np.append(y_coords, a[2]);

	y_coords.sort(key=left_coords);
	left_coords.sort();
	right_coords.sort();


	#Use the intervals themselves as testing points

	intervals = np.array(left_coords, right_coords);
	intervals = np.transpose(intervals);
	subintervals = np.sort(np.append(left_coords, right_coords));

	matrix = [];

	for a in intervals:
		matrix_line = [];
		for p in subintervals:
			if (a[0] <= p or p <= a[1]):
				matrix_line.append(1);
			else:
				matrix_line.append(0);
		matrix.append(matrix_line);

	matrix = np.transpose(np.array(matrix));

	#now that we know the intervals, integrate each circle's respective function over the interval.

	

	#End goal of this particular function is to create a list of intervals for integration
	#that will be tagged with "0, 1, 2", "2, 4", or whatever circles are intersecting them
	#After that, it will use the integration of a circle to create an analytic solution for
	#the intersection of all of those circles.

	#probably to be abandoned

