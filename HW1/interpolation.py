import numpy as np

def sorting(l1 , l2):
	l1 = np.array(l1)
	l2 = np.array(l2)
	idx = np.argsort(l1)
	return l1[idx] , l2[idx]
	
def piecewise_linear(x , y):
	'''
	This is my implementation for a piecewise linear interpolation
	x and y should both be lists / arrays that include the data.
	returns a function f(x)
	'''
	
	### first we want so sort the data in x ###
	x , y = sorting(x , y)
	
	def f(z):
		'''
		resulting function from the interpolation
		z is the point at which you would like to interpolate
		'''
		
		if z > x[-1]: ##Greater than largest x value
			print ("Warning , your data point is outside the range of available data")
			m = (y[-1] - y[-2]) / (x[-1] - x[-2])
			return y[-1] + m * (z - x[-1])
			
			
		if z < x[0]: ##Less than smallest x value
			print ("Warning , your data point is outside the range of available data")
			m = (y[0] - y[1]) / (x[0] - x[1])
			return y[0] + m * (z - x[0])			
			
		for i in range(len(x)):
			
			if z > x[i] and z < x[i+1]:
				m = (y[i] - y[i + 1]) / (x[i] - x[i + 1])
				return y[i] + m * (z - x[i])
			
				
		
	return f
	
