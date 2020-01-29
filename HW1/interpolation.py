import numpy as np
import matplotlib.pyplot as plt

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
			if z == x[i]:
				return y[i]
			if z > x[i] and z < x[i+1]:
				m = (y[i] - y[i + 1]) / (x[i] - x[i + 1])
				return y[i] + m * (z - x[i])
		
	return f
	
def h(x , i):
	return x[i + 1] - x[i]
	

def natural_cubic(x , y):
	
	### first we want so sort the data in x ###
	x , y = sorting(x , y)
	
	### Hard part is to caalculate the coefficients for all of the splines.
	### I plan to use the numpy.linalg package to solve the matrices
	
	B = [] ###matrix of b_i
	n = len(x) - 1
	
	i = 1 ###index for the b_i that we are calculating

		
	###B is set up now, so we work on A now. Our matrix problem is A * x = B
	A = []
	
	i = 1 ##row numb.
	while i < n:
	
		bi = 6 * (h(y , i) / h(x , i) - h(y , i-1) / h(x , i-1))
		B.append(bi)
		
		row = []
		k = 1 ###column numb
		while k < n:
			if k == i:
				row.append(2 * h(x , i-1) + 2 * h(x , i))
				
			elif i == k - 1:
				row.append(h(x , i))
			elif i == k + 1:
				row.append(h(x , i - 1))
			else:
				row.append(0)
			k += 1
				
		A.append(row)
		i += 1
	
	M = np.linalg.solve(A , B)
	M = np.insert(M , 0 , 0)
	M = np.append(M , 0)

	def f(z):
		#print (len(B) , len(A) , len(A[0]) , len(X))
		if z > x[-1]:
			print ("Warning , your data point is outside the range of available data")
			i = len(x) - 1
			t = z - x[i]
			a = (M[i + 1] - M[i]) / (6 * h(x , i))
			b = M[i] / 2
			c = (y[i + 1] - y[i]) / (h(x , i)) - (M[i+1] + 2 * M[i]) * h(x , i) / 6.0
			d = y[i]
			return a * t ** 3 + b * t ** 2 + c * t + d
		
		if z < x[0]:
			print ("Warning , your data point is outside the range of available data")
			i = 0
			t = z - x[i]
			a = (M[i + 1] - M[i]) / (6 * h(x , i))
			b = M[i] / 2
			c = (y[i + 1] - y[i]) / (h(x , i)) -  (M[i+1] + 2 * M[i]) * h(x , i) / 6.0
			d = y[i]
			return a * t ** 3 + b * t ** 2 + c * t + d
		
		for i in range(len(x)):
			if z  >= x[i] and z < x[i+1]:
				t = z - x[i]
				a = (M[i + 1] - M[i]) / (6 * h(x , i))
				b = M[i] / 2
				c = (y[i + 1] - y[i]) / (h(x , i)) - (M[i+1] + 2 * M[i]) * h(x , i) / 6.0
				d = y[i]

			
				
		return a * t ** 3 + b * t ** 2 + c * t + d
	return f
	
