import unittest
import numpy as np

		
		
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
	
