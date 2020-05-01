import numpy as np
import MCMC
import scipy.stats as stats
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const
##This sctipt will compute the mass and radius of our planet
#Uses the results of our lightcurve analysis.



def read_rv():
	fname = "RV_data.txt"
	f = open(fname)
	time = []
	rv = []
	
	
	for i in f.readlines():
		if i[0] == "#":
			continue
			
		time.append(float(i.split()[0]))
		rv.append(float(i.split()[1]))
		
	
	
	return time , rv
	
def Prior(X):
	
	A = stats.norm.pdf(X[0] , 3.55 , .01)
	
	return 1
	
def lsq(X , D):
	t = D[0]
	rv = D[1]
	
	ls = 0
	for i in range(len(t)):
		ls += ( X[2] * np.sin(t[i] * 2 * np.pi / X[0] + X[1]) - rv[i]) ** 2
	
	return -1 * ls * Prior(X)
	
def Q(X):

	

	cov = [ [0.08 , 0  , 0] , [0, 1e-1 ,0 ] , [0 , 0 , 10 ] ]
	mean = [0 , 0 , 0]
	Y = X + np.random.multivariate_normal(mean , cov)
	
	
	while min(Y) < 0:
		Y = X + np.random.multivariate_normal(mean , cov)
		
	return Y
	
def plot_fit(t , rv , X):
	
	fx = []
	fy = []
	
	for i in np.linspace(t[0] , t[-1] , 10000):
		fx.append(i)
		fy.append(X[2] * np.sin(i * 2 * np.pi / X[0] + X[1]))
	#print (fy)
	plt.plot(fx , fy)
	plt.scatter(t , rv)
	plt.show()

		
def measure():

	tr = 2454955.788373
	
	phi = tr * 2 * np.pi / (3.55)
	
	t , rvel = read_rv()
	

	DI = 0.005
	Rstar = 1.79
	
	Rp = np.sqrt(Rstar ** 2 * DI)
	
	print (Rp)
	
	###Kepler's Third Law
	
	#a^3 / T^2 = G(M + m) / (4pi^2)
	# pplanet = known
	#a = G M / v^2 --> a = G M / (pplanet ** 2 / m **2 ) = G * M * m ** 2 / (Pplanet ** 2)
	
	# 4 pi^2 ( GMm^2/(Pp^2) )^3 = G (M + m )
	
	
	x0 = [3.55 , 4.34505429e+06 , 2.50998989e+02]

	opt = MCMC.MCMC(lsq , [t , rvel] )
	opt.M_H(Q , 5000 , x0 )
	print (opt.MAP , opt.MAPL)
	plot_fit(t , rvel , opt.MAP)
	
	vp = opt.MAP[2] * u.m / u.s
	
	Mstar = 1.35 * u.M_sun
	
	pstar = Mstar * vp
	p_planet = pstar
	
	
	T = opt.MAP[0] * u.day
	m_planet = T ** 2 * (p_planet ** 6) / (4 * np.pi ** 2 * (const.G * Mstar) ** 2)
	m_planet = m_planet ** (1.0 /6 )
	print (m_planet.to(u.M_sun))

	
	
measure()

