import unittest
import matplotlib.pyplot as plt
import astropy.units as u

def test_ode(t , y):
	
	##solves y'(t) = 2 * t --> y = t ^ 2
	
	return [2 * t]

def test_coupled(t , y):

	K = 100
	y_0 = 1.2 * y[0] * (K - y[0]) / K
	y_1 =  y[1] * (K - y[1]) / K - 0.5 * y[0]
	return [y_0 , y_1]
	
class solve_ode:
	
	def __init__(self , f , t0 , y0 , h = 1e-3):
		self.f = f
		self.t0 = t0
		self.y0 = y0
		self.h = h
		## f is a function of y and t -- f(t , y)
	
	def Forward_Euler(self, t_end):
		
		'''
		solves our ode using a forwad euler method
		t_end should be the final time
		returns two lists, t and y
		t contains the time of each output
		y is an array containing the valuables to all of our variables
		'''
		
		
		t = []
		y = []
		
		t.append(self.t0)
		y.append(self.y0)
		
		while t[-1] < t_end:
			
			fi = self.f(t[-1] , y[-1])
			t.append(t[-1] + self.h)
			ny = []
			for i in range(len(y[-1])):
				ny.append( y[-1][i] + fi[i] * self.h)
			y.append(ny)
			
		return t , y

	def Heun(self , t_end):
		
		t = []
		y = []
		
		t.append(self.t0)
		y.append(self.y0)
		
		while t[-1] < t_end:
			
			fi = self.f(t[-1] , y[-1])
			t.append(t[-1] + self.h)
			predictor = []
			for i in range(len(y[-1])):
				predictor.append( y[-1][i] + fi[i] * self.h)
			
			###predictor is the result of the predictor step
			
			fi1 = self.f(t[-1] , predictor)
			ny = []
			for i in range(len(y[-1])):
				ny.append( y[-1][i] + 0.5 * self.h * (fi1[i] + fi[i]))
			y.append(ny)
			
		return t , y
		
def test():
	A = solve_ode(test_coupled , 0 , [5 , 10] )
	t , y = A.Heun(10)
	
	py = []
	for i in y:
		py.append(i[0])
		
	plt.plot(t , py)
	plt.show()
	
	
test()
