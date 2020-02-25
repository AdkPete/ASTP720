import unittest
import matplotlib.pyplot as plt
import astropy.units as u

	
def add_lists(y1 , y2):

	'''
	Takes in two lists of the same length, y1 and y2
	creates a new list y3 such that y3[i] = y1[i] + y2[i]
	returns y3
	'''
	y3 = []
	for i in range(len(y1)):
		y3.append(y1[i] + y2[i])
	return y3
	
def list_scale(y1 , scale):

	'''
	takes in an list y1 and a float called scale
	creates a new alist y2 such that y2[i] = y1[i] * scale
	returns y2
	'''
	
	y2 = []
	for i in range(len(y1)):
		y2.append(y1[i] * scale)
	return y2
	
class solve_ode:
	
	'''
	This class will solve an ode
	includes three different ode solvers, which are Forward Euler, Heuns method, and an RK4 scheme
	'''
	
	def __init__(self , f , t0 , y0 , h = 1e-3 , n = None , t_end = None):
		'''
		f is a function such that y'(t) = f(t , y)
		in problems with more than one dimension, f(t , y) should return an list with the same shape as y
		t0 is the time for the initial condition
		y0 is the initial condition, should be an list or float for 1D problems, should be an list for higher dimensions
		h is the step size, which defaults to 1e-3
		n is the number of steps you want to take. Setting n will overwrite h
		t_end is the end time. Can be set here or in the actual sovlers
		'''
		
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
			
		###The following will make it such that you can solve problems with float inputs instead of lists
		
		if type(y0) == type(1.01) or type(y0) == type(1):
			self.y0 = [self.y0]
			
		if type(self.f(t0 , y0)) == type(1.01) or type(self.f(t0 , y0)) == type(1):
			orig_func = self.f
			def new_func(t , y):
				return [ orig_func(t , y) ]
			self.f = new_func
		
		
	def set_h(self):
		'''
		sets h if we have specified a number of steps
		'''
		
		if self.use_n: ###here we will calculate the step size based on the number of steps
			self.h = (self.t_end - self.t0) / self.n
			
	def clean_output(self , y):
		'''
		Does nothing if there is more than one variable
		in the one variable case, it turns each element of y from an list to a float
		'''
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
		y is an list containing the values to all of our variables
		also sets self.fe_t = t and self.fe_y = y
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
		y is an list containing the values to all of our variables
		also sets self.heun_t = t and self.heun_y = y
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
			
			f_i = self.f(t[-1] , y[-1])
			
			###Starts with predictor step
			
			predictor = add_lists(list_scale(f_i, self.h) , y[-1])
			
			###Now for the corrector step
			
			corr = add_lists( f_i , self.f(t[-1] + self.h , predictor))
			
			corrected = add_lists(y[-1] , list_scale(corr , 0.5 * self.h))
			
			
			y.append(corrected)
			t.append(t[-1] + self.h)
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
		y is an list containing the values to all of our variables
		also sets self.rk4_t = t and self.rk4_y = y
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
			
			k1 = list_scale(self.f(t[-1] , y[-1]) , self.h)
			k2 = list_scale(self.f(t[-1] + self.h / 2.0 , add_lists(y[-1] , list_scale(k1 , 0.5))) , self.h)
			k3 = list_scale(self.f(t[-1] + self.h / 2.0 , add_lists(y[-1] , list_scale(k2 , 0.5))) , self.h)
			k4 = list_scale(self.f(t[-1] + self.h, add_lists(y[-1] , k3)) , self.h)

			rk = add_lists(k1 , list_scale(k2 , 2))
			rk = add_lists(rk , list_scale(k3 , 2))
			rk = add_lists(rk , k4)
			rk = list_scale(rk , (1 / 6.0))
			
			t.append(t[-1] + self.h)
			y.append(add_lists(y[-1] , rk))
			
			
		y = self.clean_output(y)
		
		self.rk4_t = t
		self.rk4_y = y
		return t , y
		
