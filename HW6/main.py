import numpy as np


sigma = 1

def in_dip(tref , tau , ti , P):
	
	n = 1
	while True:
		time = ti - n * p
		if time < tref:
			return False
		elif time > tref and time < tref + tau:
			return True
		n += 1
	
		
def s(tref , tau ,  ti , DI , P):
	
	'''
	
	This is the periodic dip function
	Takes in a variety of parameters
	Returns a change in intensity
	'''
	
	
	if in_dip(tref , tau , ti , P):
		return DI
	else:
		return 0

def L(t , I , tref , tau , DI , P):
	'''
	Returns the value for our likelihood function
	ti and I are arrays containing our light curve data
	'''
	
	for i in range(len(I)):
		if i == 0:
			A = np.exp(-0.5 * (1 / sigma ** 2) * (I[i] - s(tref , tau , t[i] , DI , P) ** 2
		else:
			A *= np.exp(-0.5 * (1 / sigma ** 2) * (I[i] - s(tref , tau , t[i] , DI , P) ** 2
		
	return 1 / np.sqrt(2 * np.pi * sigma * sigma) * A
	

	
	
	
