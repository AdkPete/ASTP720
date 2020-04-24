import numpy as np
import MCMC
from scipy.stats import multivariate_normal
import matplotlib.pyplot as plt

def poisson( lam , k):
	return lam ** k * np.exp(-1 * lam) / (np.math.factorial(k))
	
def Q0(X):
	cov = [ [ .5 , 0 ] , [0 , .5 ] ] 
	mean = [0 , 0 ]
	Y = X + np.random.multivariate_normal(mean , cov)
	return Y
	
def P0(X , data):
	sig = [ [1.5 , 0 ], [0 , 3.4 ]]
	mu = [ 3 , 10 ]
	rv = multivariate_normal( mu , sig )
	return rv.pdf(X)
	

def Q1(X):
	A = np.random.rand()
	if X[0] < 1:
		return X + [1]
	
	if A >= 0.5:
		return X + [1]
	else:
		return X - [1]
		
def P1(X , data):

	offset = 1e4
	L = 0
	for i in range(len(data)):
		A = poisson(X[0] , data[i])
		if A == 0:
			A = 10 ** -50
		L += np.log10( A )
		
	return L + offset
		
	
def mkp0(X):

	x0 = []
	x1 = []
	for i in X:
		x0.append(i[0])
		x1.append(i[1])
		
	
	
	plt.subplot(2 , 1 , 1)
	
	plt.hist(x0)
	
	plt.subplot(2 , 1 , 2)
	
	plt.hist(x1)
	
	plt.show()
	
def mkp1(X):
	x0 = []
	for i in X:
		
		x0.append(i[0])
	
	plt.hist(x0)
	plt.show()
	
def T0():
	data = 2
	T = MCMC.MCMC( P0 , data)


	X, MP , ML = T.M_H(Q0 , 10000 , np.array(  [ 10 , 3 ] ) )
	mkp0(X)
	
def T1():
	
	data = np.random.poisson(3 , 1000)
	T = MCMC.MCMC(P1 , data)
	
	X , MP , ML = T.M_H(Q1 , 10000 , np.array( [ 5 ] ) )
	mkp1(X)
	print (MP , ML)
	
#T0()
T1()
