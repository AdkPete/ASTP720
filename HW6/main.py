import numpy as np
import MCMC

sigma = 1

def in_dip(tref , tau , ti , P):
	
	n = 1
	while True:
		time = ti - n * P
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
			A = -0.5 * (1 / sigma ** 2) * np.exp(I[i] - s(tref , tau , t[i] , DI , P)) ** 2
		else:
			A *= -0.5 * (1 / sigma ** 2) * np.exp(I[i] - s(tref , tau , t[i] , DI , P)) ** 2
		
	return 1 / np.sqrt(2 * np.pi * sigma * sigma) * A
	
def P(X , Data):
	
	t = Data[0]
	I = Data[1]
	
	tref = X[0]
	tau = X[1]
	DI = X[2]
	P = X[3]
	
	return L(t , I , tref , tau , DI , P)
	
def Q(X):
	Y = []
	for i in range(len(X)):
		Y.append(X[i] + np.random.normal(0 , 1 , 1)[0])
	return Y
	
def read_lightcurve(fname):
	f = open(fname)
	
	time = []
	flux = []
	for i in f.readlines():
		if i[0] == "#":
			continue
			
		time.append(float(i.split()[0]))
		flux.append(float(i.split()[1]))
		
	return time , flux
	
t , f = read_lightcurve("lightcurve_data.txt")
data = [t , f]
A = MCMC.MCMC(P , data)
A.M_H(Q , 10 , [t[0] , 1 , 1 , 30])
	

		
	
