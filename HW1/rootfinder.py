import sys
import numpy as np


def test(x):
	'''
	takes in a value x
	returns a polynomial evaluated at x
	'''
	
	return x ** 3 - 5 * x ** 2 + 17 * x - np.sqrt(133)

def testprime(x):
	'''
	takes in a point x
	returns the derivative of our test function at that point
	'''
	return 3 * x ** 2 - 10 * x + 17

def bisection(f , xlow , xhigh , threshhold =  1e-5 , iterations = False):
	##Implementation of the bisection root finding algorithm
		
	'''
	
	f should be the function you want to find a root of
	xlow and xhigh are the two starting boundaries for the bisection algorithm.
	There shoule be one root between the two of these
	threshhold sets the accuracy threshold
	iterations determines wether or not the code will output the number of iterations
	
	'''
	
	if f(xlow) * f(xhigh) > 0: 
		###Bisection will not work if there is no root inside of your starting interval
		print ("Error , even number of roots in interval, or 0 roots in interval")
		sys.exit()
	
	niter = 0
	
	while (xhigh - xlow) > threshhold:
		xmid = (xhigh + xlow) / 2.0
		
		if f(xmid) * f(xhigh) < 0:  ###root is in upper half of interval
			xlow = xmid
			
		else: 						###root is in lower half of interval
			xhigh = xmid
		niter += 1
			
	if iterations:
		return (xhigh + xlow) / 2.0 , niter
			
	return (xhigh + xlow) / 2.0
		
		
def Newton(f , fprime , x , threshhold = 1e-5 , iterations = False):
	####Implementation of the newton method root finding algorithm
	
	'''
	f should be the function you want to find a root of
	fprime is the derivative of f
	x is your starting guess
	threshhold sets the accuracy threshold
	iterations determines wether or not the code will output the number of iterations
	'''
	
	niter = 0
	
	while abs(f(x)) > threshhold:
		niter += 1
		x = x - f(x) / fprime(x) ###Determines the next x value
	
	if iterations:
		return x , niter
	return x
	
def Secant(f  , x0 , x1 , threshhold = 1e-5 , iterations = False):
	####Implementation of the secant method root finding algorithm

	'''
	f should be the function you want to find a root of
	x0 and x1 are your starting guessess
	threshhold sets the accuracy threshold
	iterations determines wether or not the code will output the number of iterations
	'''
	
	niter = 0
	
	while abs(f(x1)) > threshhold:
		niter += 1
		x = x1 - (f(x1) * (x1 - x0)) / (f(x1) - f(x0)) ###Determines the next x value
		x0 = x1
		x1 = x
		
	if iterations:
		return x1 , niter
		
	return x1
		
def run_test():
	'''
	Tests all three root finding algorithms
	prints out end results and the number of iterations requires
	'''
	numerical , niter = bisection(test , -1 , 1 , 1e-10 , True)
	print (numerical , niter)
	print (test(numerical))
	
	
	numerical , niter = Newton(test , testprime , 2 , 1e-10 , iterations = True)
	print (numerical , niter)
	print (test(numerical))
	
	numerical , niter = Secant(test , -2 , 2 ,  1e-10 , iterations = True)
	print (numerical , niter)
	print (test(numerical))

