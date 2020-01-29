import numpy as np
import matplotlib.pyplot as plt

def sorting(l1 , l2):
	'''
	Takes in two lists or np arrays
	Sorts them in order of the first list.
	Returns 2 sorted np arrays
	'''
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
	##takes in a list and position, returns the spacing between adjacent elements
	return x[i + 1] - x[i]
	

def natural_cubic(x , y):
	'''
	This is my implementation for natural cubic spline interpolation
	Takes in two lists containing all of the data
	Returns a function, f(x) which will interpolate at a point x
	'''
	
	### first we want so sort the data in x ###
	x , y = sorting(x , y)
	
	###First we compute all of the matrices
	
	B = [] ###matrix of b_i
	n = len(x) - 1
	

	A = [] ### this will be our H_ij
	
	i = 1 ##row number
	while i < n:
	
		##Note that the g_i that appear in the b_i can be calculated using h(y , i)
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
	
	###M contains all of our second derivatives

	def f(z):
		#print (len(B) , len(A) , len(A[0]) , len(X))
		if z > x[-1]:
			###We are extrapolating here
			print ("Warning , your data point is outside the range of available data")
			i = len(x) - 1
			
			### It is useful to write the splines in terms of t = x - x_i
			t = z - x[i]
			
			a = (M[i + 1] - M[i]) / (6 * h(x , i))
			b = M[i] / 2
			c = (y[i + 1] - y[i]) / (h(x , i)) - (M[i+1] + 2 * M[i]) * h(x , i) / 6.0
			d = y[i]
			
			
			return a * t ** 3 + b * t ** 2 + c * t + d
		
		if z < x[0]:
			###We are extrapolating here
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
	
