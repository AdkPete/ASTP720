import astropy , sys
import numpy as np


def test(x):
	return x ** 3 - 5 * x ** 2 + 17 * x - np.sqrt(133)

def testprime(x):
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
	
	if f(xlow) * f(xhigh) > 0: ###Bisection will not work if there is no root inside of your starting interval
		print ("Error , even number of roots in interval, or 0 roots in interval")
		sys.exit()
	
	niter = 0
	
	while (xhigh - xlow) > threshhold:
		xmid = (xhigh + xlow) / 2.0
		
		if f(xmid) * f(xhigh) < 0:
			xlow = xmid
			
		else:
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
		x = x - f(x) / fprime(x)
	
	if iterations:
		return x , niter
	return x
		
def run_test():
	numerical , niter = bisection(test , -1 , 1 , 1e-10 , True)
	print (numerical , niter)
	print (test(numerical))
	
	
	numerical , niter = Newton(test , testprime , 20 , 1e-10 , iterations = True)
	print (numerical , niter)
	print (test(numerical))
	
run_test()
	
	

