import unittest
import numpy as np


class TestCalc(unittest.TestCase):
	
	#runs a number of test cases
	
	def test_centered_diff(self):
		self.assertAlmostEqual( 5 , centered_diff(test_func_1 ,1e-4)(5))
		
	def test_rectangle(self):
		###note the low accuracy restriction due to the inefficient nature of this method
		self.assertEqual( round(right_rect_integrate(test_func_1 , 0 , 5 , 1e-5) , 3) , 112.5)
		
	def test_trap(self):
		
		self.assertAlmostEqual( trap_rule(test_func_1 , 0 , 5 , 1e-4) , 112.5)
		self.assertAlmostEqual( trap_rule(test_func_2 , 0 , np.pi * 2 , 1e-4) , 0)
		
	def test_midpoint(self):
		self.assertAlmostEqual( midpoint_rule(test_func_1 , 0 , 5 , 1e-4) , 112.5)
		self.assertAlmostEqual( midpoint_rule(test_func_2 , 0 , np.pi * 2 , 1e-4) , 0)
		
	def test_simpson(self):
		self.assertAlmostEqual( simpson(test_func_1 , 0 , 5 , 1e-4) , 112.5)
		self.assertAlmostEqual( simpson(test_func_2 , 0 , np.pi * 2 , 1e-4) , 0)
		
		
def test_func_1(x):
	'''
	This is a test function
	'''
	
	return 5 * x + 10
	
def test_func_2(x):
	'''
	This is a second test function
	'''
	
	return np.sin(x)


def centered_diff(f , h):
	
	'''
	this function will return the derivative of the function f at a point x
	f should be a function
	x and h both are numbers
	returns a function for the derivative
	'''
	
	return lambda x: (f(x + h) - f(x - h)) / (2 * h)
	
def right_rect_integrate(f , a , b , h):

	'''
	this function will return the integral of f from a to b
	uses a rectangle method
	f should be a function
	a and b are the limits of integration
	h is the width of the rectangles
	returns the value of the integral
	'''
	
	
	x = a
	integral = 0
	while x < b - h:
		integral += f(x) * h
		x += h
	
	###handles any remaining bit of the region
	if x < b:
		nh = b - x
		integral += f(x) * nh

	return integral
	
def trap_rule(f , a , b , h):

	'''
	This function will integrate using the trap rule
	a and b are the limits of integration
	f is the function that you would like to integrate
	h is the width of the trapezoids
	returns the value of the integral
	'''
	
	x = a
	integral = 0
	while x < b - h:
		integral += (f(x) + f(x + h)) * h / 2
		x += h
		
	if x < b:
		nh = b - x
		integral += (f(x) + f(x + nh)) * nh / 2
		
	return integral
	
def midpoint_rule(f , a , b , h):
	'''
	This function will integrate using the midpoint rule
	a and b are the limits of integration
	f is the function that you would like to integrate
	h is the width of each rectangle
	returns the value of the integral
	'''
	
	x = a
	integral = 0
	while x < b - h:
		integral += f((2 * x + h) / 2) * h
		x += h
		
	if x < b:
		nh = b - x
		integral += f((2 * x + nh) / 2) * nh
		
	return integral
	
def simpson(f , a , b , h):
	'''
	This function will integrate using the Simpson rule
	a and b are the limits of integration
	f is the function that you would like to integrate
	h is the width of each region
	returns the value of the integral
	'''
	x = a
	integral = 0
	while x < b - h:
		
		integral += (h / 6) * (f(x) + f(x + h) + 4 * f((2 * x + h) / 2))
		x += h
		
	if x < b:
		nh = b - x
		integral += (nh / 6) * (f(x) + f(x + nh) + 4 * f((2 * x + nh) / 2))
		
	return integral
	
if __name__ == "__main__":
	'''
	If we run this script, a series of tests will run to check
	all of the important functions in this file
	'''
	unittest.main()
