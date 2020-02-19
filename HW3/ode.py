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
	
def add_arrays(y1 , y2):
	y3 = []
	for i in range(len(y1)):
		y3.append(y1[i] + y2[i])
	return y3
	
def array_scale(y1 , scale):
	y2 = []
	for i in range(len(y1)):
		y2.append(y1[i] * scale)
	return y2
	
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
		y is an array containing the values to all of our variables
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
		
		'''
		solves our ode using Heun's method
		t_end should be the final time
		returns two lists, t and y
		t contains the time of each output
		y is an array containing the values to all of our variables
		'''
		
		
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
		
	def RK4(self , t_end):
		
		'''
		solves our ode using an RK4 scheme
		t_end should be the final time
		returns two lists, t and y
		t contains the time of each output
		y is an array containing the values to all of our variables
		'''
		
		t = []
		y = []
		
		t.append(self.t0)
		y.append(self.y0)
		
		while t[-1] < t_end:
			
			k1 = self.f(t[-1] , y[-1])
			y2 = add_arrays(y[-1] , array_scale(k1 , 0.5))
			
			k2 = self.f(t[-1] + 0.5 * self.h , y2)
			
			y3 = add_arrays(y[-1] , array_scale(k2 , 0.5))
			k3 = self.f(t[-1] + 0.5 * self.h , y3)
			y4 = add_arrays(y[-1] , k3)
			k4 = self.f(t[-1] + self.h , y4)
			
			t.append(t[-1] + self.h)
			k12 = add_arrays(k1 , k2)
			k123 = add_arrays(k12 , k3)
			k1234 = add_arrays(k123 , k4)
			rk = array_scale(k1234 , self.h * 1.0 / 6)
			
			y.append(add_arrays(y[-1] , rk))
			
		return t , y
		
		
def test():
	A = solve_ode(test_coupled , 0 , [5 , 10] )
	t , y = A.RK4(10)
	py = []
	for i in y:
		py.append(i[0])
		
	plt.plot(t , py)
	plt.show()
	
	
if __name__ == "__main__":
	test()
