import unittest
import matplotlib.pyplot as plt
import astropy.units as u

	
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
	
	def __init__(self , f , t0 , y0 , h = 1e-3 , n = None , t_end = None):
		self.f = f
		self.t0 = t0
		self.y0 = y0
		self.h = h
		self.use_n = False
		self.n = n
		
		self.fe_t = None
		self.fe_y = None
		self.heun_t = None
		self.heun_y = None
		self.rk4_t = None
		self.rk4_y = None
		
		if t_end != None:
			self.t_end = t_end
			
		if n != None:
			self.use_n = True
			
		###The following will make it such that you can solve problems with float inputs instead of arrays
		
		if type(y0) == type(1.01) or type(y0) == type(1):
			self.y0 = [self.y0]
			
		if type(self.f(t0 , y0)) == type(1.01) or type(self.f(t0 , y0)) == type(1):
			orig_func = self.f
			def new_func(t , y):
				return [ orig_func(t , y) ]
			self.f = new_func
		
		
	def set_h(self):
		if self.use_n: ###here we will calculate the step size based on the number of steps
			self.h = (self.t_end - self.t0) / self.n
			
	def clean_output(self , y):
		if len(y[0]) == 1:
			new_y = []
			for i in y:
				new_y.append(i[0])
			return new_y
		return y
		
	
	def Forward_Euler(self, t_end = None):
		
		'''
		solves our ode using a forwad euler method
		t_end should be the final time
		returns two lists, t and y
		t contains the time of each output
		y is an array containing the values to all of our variables
		'''
		
		if t_end == None:
			t_end = self.t_end
		else:
			self.t_end = t_end
		
		self.set_h()
		
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
			
		y = self.clean_output(y)
		self.fe_t = t
		self.fe_y = y
		return t , y

	def Heun(self , t_end = None):
		
		'''
		solves our ode using Heun's method
		t_end should be the final time
		returns two lists, t and y
		t contains the time of each output
		y is an array containing the values to all of our variables
		'''
		
		if t_end == None:
			t_end = self.t_end
		else:
			self.t_end = t_end
		
		self.set_h()
		
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
			
		y = self.clean_output(y)
		
		self.heun_t = t
		self.heun_y = y
		return t , y
		
	def RK4(self , t_end = None):
		
		'''
		solves our ode using an RK4 scheme
		t_end should be the final time
		returns two lists, t and y
		t contains the time of each output
		y is an array containing the values to all of our variables
		'''
		
		if t_end == None:
			t_end = self.t_end
		else:
			self.t_end = t_end
		
		self.set_h()
		
		t = []
		y = []
		
		t.append(self.t0)
		y.append(self.y0)
		
		while t[-1] < t_end:
			
			k1 = array_scale(self.f(t[-1] , y[-1]) , self.h)
			k2 = array_scale(self.f(t[-1] + self.h / 2.0 , add_arrays(y[-1] , array_scale(k1 , 0.5))) , self.h)
			k3 = array_scale(self.f(t[-1] + self.h / 2.0 , add_arrays(y[-1] , array_scale(k2 , 0.5))) , self.h)
			k4 = array_scale(self.f(t[-1] + self.h, add_arrays(y[-1] , k3)) , self.h)

			rk = add_arrays(k1 , array_scale(k2 , 2))
			rk = add_arrays(rk , array_scale(k3 , 2))
			rk = add_arrays(rk , k4)
			rk = array_scale(rk , (1 / 6.0))
			
			t.append(t[-1] + self.h)
			y.append(add_arrays(y[-1] , rk))
			
			
		y = self.clean_output(y)
		
		self.rk4_t = t
		self.rk4_y = y
		return t , y
		
