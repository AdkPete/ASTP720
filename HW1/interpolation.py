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
			if z == x[i]:
				return y[i]
			if z > x[i] and z < x[i+1]:
				m = (y[i] - y[i + 1]) / (x[i] - x[i + 1])
				return y[i] + m * (z - x[i])
		
	return f

def natural_cubic(x , y):
	
	### first we want so sort the data in x ###
	x , y = sorting(x , y)
	
	### Hard part is to caalculate the coefficients for all of the splines.
	### I plan to use the numpy.linalg package to solve the matrices
	
	B = [] ###matrix of b_i
	n = len(x) - 1
	
	i = 1 ###index for the b_i that we are calculating
	
	while i < n:
		hi = x[i + 1] - x[i]
		him1 = x[i] - x[i - 1]
		gi = y[i + 1] - y[i]
		gim1 = y[i] - y[i-1]
		
		bi = 6 * (gi / hi - gim1 / him1)
		B.append(bi)
		i += 1
		
	###B is set up now, so we work on A now. Our matrix problem is A * x = B
	A = []
	
	i = 1 ##row numb.
	while i < n:
		row = []
		k = 1 ###column numb
		while k < n:
			if k == i:
				him1 = him1 = x[i] - x[i - 1]
				hi = x[i + 1] - x[i]
				row.append(2 * him1 + hi)
				
			elif i == k - 1:
				hi = x[i + 1] - x[i]
				row.append(hi)
			elif i == k + 1:
				him1 = x[i] - x[i - 1]
				row.append(him1)
			else:
				row.append(0)
			k += 1
				
		A.append(row)
		i += 1
	
	X = np.linalg.solve(A , B)
	
	def f(z):
		print (len(B) , len(A) , len(A[0]) , len(X))
		return z
		
	return f
	
def g(x):
	return x ** 2
	
def test():
	x = range(-10 , 11)
	y = []
	
	for i in range(len(x)):
		y.append(g(x[i]))
		
	f = natural_cubic(x , y)
	
	print(f(4.7))
	
test()
