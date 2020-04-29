import numpy as np
import MCMC
import matplotlib.pyplot as plt
import sys

#####TODO######

'''
4. Solve with MCMC
5. Solve with Genetic Algorithm
'''

def fold(time , flux , period , tref , tau):
	
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

	cov = [ [ .01 , 0 , 0 , 0 ] , [0 , .01 , 0 , 0] , [0 , 0 , .00001 , 0 ] , [0 , 0 , 0 , 0.0001 ] ] 
	mean = [0 , 0 , 0 , 0]
	Y = X + np.random.multivariate_normal(mean , cov)
	
	
	while min(Y) < 0:
		Y = X + np.random.multivariate_normal(mean , cov)
	return Y
	
def s(period, DI , tau , t):
	
	start = period / 2.0
	if t >= start and t <= tau + start:
		
		return 1 - DI
	
	return 1
	
def mkplot(X , D):
	nt , nf = fold(t , f , X[0] , X[1] , X[3])
	
	
	tt = []
	tf = []
	for i in np.linspace(0 , max(nt) , 300):
		tt.append(i)
		tf.append(s(X[0] , X[2] , X[3] , i ))
		
	plt.plot(tt , tf , color = "green" )
	plt.scatter(nt , nf)
	plt.show()
	

def L(X , D):
	time = D[0]
	flux = D[1]
	
	P = 3.55
	time , flux = fold(time , flux , P , X[1] , X[3])
	

	period = P
	tref = X[1]
	DI = X[2]
	tau = X[3]
	
	sigma = 0.002
	
	Lhood = 0
	
	for i in range(len(time)):
	
		A = 1 / np.sqrt(2 * np.pi * sigma)
		B = -0.5 * (1 / (sigma ** 2)) * (flux[i] - s(period, DI , tau , time[i])) ** 2
		if A * np.exp(B) == 0 or tau < 0.1 or DI < .004 or period > 4 or period < 3 or DI > .01 or tau > 2:
			Lhood += -30
		else:
			Lhood += np.log10(A * np.exp(B))
	
	return Lhood

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

rx = [ 3.55 , 2.2 , .007 , 0.15 ]

x0 = [3 , 1 , 0.017272 , 0.32]
bl = L(x0 , [t , f] )


for i in range(100):
	print (x0 , bl)
	opt = MCMC.MCMC(L , [t , f])
	opt.M_H(Q , 1000 , x0)
	if opt.MAPL > bl:
		bl = opt.MAPL
		x0 = opt.MAP


print (opt.MAP  , opt.MAPL )

mkplot( opt.MAP , [t , f] )
