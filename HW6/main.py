import numpy as np
import MCMC
import matplotlib.pyplot as plt
import sys , time
from scipy.stats import norm

###This script will estimate the paremeters for our light curve
###Comes up with the period, tref

def avg_lcurve(ft , fflux):
	
	'''
	Takes in a folded light curve
	Averages together nearby points
	Returns new averaged light curve values
	'''
	avg_t = []
	avg_f = []
	
	for i in range(0 , len(ft) , 5):
		avg_t.append(np.mean( ft[i:i+5] ) )
		avg_f.append(np.mean( fflux[i:i+5] ) )
		

		
	return avg_t , avg_f
	

def fold(time , flux , period , tref , tau):
	
	'''
	Code to fold a lightcurve
	time and flux should be arrays of light curve values
	the period, tref and tau are all parameters for the transit
	Returns two arrays, one for the folded time and one for the folded flux
	'''
	
	t_offset = time[0]
	for i in range(len(time)):
		time[i] = time[i] - t_offset
	
	
	
	fold_time = []
	fold_flux = []
	start = tref - period / 2.0 + tau / 2.0
	n = 0
	
	for i in range(len(time)):
		if time[i] < start:
			continue
			
		#print (time[i] - 0 * period - start)

		if (time[i] - n * period - start) < period and (time[i] - n * period - start) > 0:
			fold_time.append(time[i] - n * period - start)
			fold_flux.append(flux[i])
		else:
			n += 1
			
	
	return fold_time , fold_flux
	
def Q(X):

	'''
	A proposal distribution for solving all 4 parameters
	Tuned to provide a sensible acceptance rate
	'''

	cov = [ [0.08 , 0 , 0 , 0 ] , [0, 1e-1 , 0 , 0 ] , [0 , 0,  1e-3 , 0 ]  , [ 0 , 0 , 0 , .1 ]] 
	mean = [0 , 0 , 0 , 0]
	Y = X + np.random.multivariate_normal(mean , cov)
	
	
	while min(Y) < 0:
		Y = X + np.random.multivariate_normal(mean , cov)
	return Y
	
def Q2(X):
	'''
	A different proposal distribution, used for testing only
	'''
	

	cov = [[1 , 0 , 0 , 0], [ 0 , 1e-50 , 0 , 0 ] , [0 , 0 , 1e-2 , 0 ]  , [ 0 , 0 , 0 ,1e-50 ]] 
	mean = [ 0 , 0 , 0 , 0]
	Y = X + np.random.multivariate_normal(mean , cov)
	
	
	while min(Y) < 0:
		Y = X + np.random.multivariate_normal(mean , cov)
	return Y
	
def s(period, DI , tau , t):

	'''
	Takes in lightcurve parameters, period, Delta I, and tau
	part of the lilelihood function
	if we are in a transit, returns 1 - DI
	otherwise, returns 1
	'''
	
	start = period / 2.0
	if t >= start and t <= tau + start:
		
		return 1 - DI
	
	return 1
	
def plot_fit(X , D):

	'''
	plots both a scatter plot of our data and our fit
	produces a plot, no returns
	'''
	
	nt , nf = fold(D[0] , D[1] , X[0] , X[1] , X[3])
	
	nt , nf = avg_lcurve(nt , nf)
	
	tt = []
	tf = []
	for i in np.linspace(0 , max(nt) , 300):
		tt.append(i)
		tf.append(s(X[0] , X[2] , X[3] , i ))
		
	plt.plot(tt , tf , color = "green" )
	plt.scatter(nt , nf)
	plt.show()
	
def Prior(X):
	
	##Prior distribution
	##Keeps parameters within reasonable ranges to avoid local maxima
	
	tref = X[1]
	DI = X[2]
	tau = X[3]
	if  tau < 0.1 or DI < .004 or X[0] > 5 or X[0] < 2 or tau > .5 or tref > 3 or tref < 1:
		##THis part of the prior asserts that there must be a transit, so setting delta_i = tau = 0 is not a valid solution
		return 1e-75
	return 1
	

'''
def L(X , D):

	
	#Likelihood function
	#Takes in our parameters X and data D, returns the value of the likelihood
	
	time = D[0]
	flux = D[1]
	
	time , flux = fold(time , flux , X[0] , X[1] , X[3])
	time , flux = avg_lcurve(time , flux)

	period = X[0]
	tref = X[1]
	DI = X[2]
	tau = X[3]
	
	sigma = 1
	
	Lhood = 1
	
	for i in range(len(time)):
	
		A = 1 / np.sqrt(2 * np.pi * sigma)
		B = -0.5 * (1 / (sigma ** 2)) * (flux[i] - s(period, DI , tau , time[i])) ** 2
		
		if A * np.exp(B) == 0:
			return 0
		else:
			Lhood *= A * np.exp(B) * 2.5
	
	return Lhood + np.log10(Prior(X))
'''

def P(X , D):
	##Similar to above, this is a different evaluation function
	##I wrote this up as an experiment
	t = D[0]
	f = D[1]
	P = X[0]
	tref = X[1]
	DI = X[2]
	tau = X[3]
	sigma = 1
	nt , nf = fold(t , f , P , tref , tau)
	
	nt , nf = avg_lcurve(nt , nf)
	A = 1.0 / np.sqrt(2 * np.pi * sigma)
	Lhood = 0.01
	for i in range(len(nf)):
		
		Lhood += (nf[i] - s(P , DI , tau , nt[i])) ** 2
		
	return (1 / Lhood) * Prior(X)
		


def read_lightcurve(fname):
	##Reads a lightcurve from a file named fname
	##Returns arrays of times and fluxes
	
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


x0 = [3.5 , 2.2 , .1 , 0.3] ##A reasonable initial guess

r = [3.546 , 2.25 , .007 , .15] ##This is my current estimate for our best solution


opt = MCMC.MCMC(P , [t , f])


X , M , ML = opt.M_H(Q2 , 5000 , r) ##Lots of iterations, so somewhat slow. Good chance to find the right answer.


print (opt.MAP , opt.MAPL)
plot_fit(opt.MAP , [t , f])

p_list = []
di_list = []
for i in X[500::]:
	
	p_list.append(i[0])
	di_list.append(i[2])
	

plt.subplot(1 , 2 , 1)
plt.hist(p_list , bins = 25)
plt.xlabel("Period")
plt.ylabel("Number of Steps")

plt.subplot(1 , 2 , 2)
plt.hist(di_list , bins = 25)
plt.xlabel("Delta I")
plt.ylabel("Number of Steps")
plt.savefig("Hist.pdf")

