import astropy , sys
import numpy as np


def test(x):
	return x ** 2 - 150

def bisection(f , xlow , xhigh , threshhold =  .0001 , iterations = False):
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
		
		
def run_test():
	correct = np.sqrt(150)
	numerical , niter = bisection(test , 10 , 20 , 1e-10 , True)
	print (numerical , niter)
	print (numerical - correct)

	
	

