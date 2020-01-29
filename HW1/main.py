import rootfinder
import interpolation
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy import constants as const

def	f(t):
	return (1 + t ** 2) ** (-0.5) - .5

def fprime(t):
	return -1 * t / ((1 + t ** 2) ** 1.5)
	
def lens(x , wlength , N0 , D , a):
	#re = const.e.esu ** 2 / (const.m_e * const.c ** 2)
	re = 2.817940328e-15 * u.m
	C = ((wlength ** 2 * N0 * re * D) /(np.pi * a * a)) * np.exp(-((x / a) ** 2))
	#print ((C).to(u.dimensionless_unscaled))
	return x * (1 + C)
	
def problem_3():
	
	thresh = 1
	th = []
	bi_iter = []
	s_iter = []
	n_iter = []
	while thresh > 1e-10:
		b_value , b_niter = rootfinder.bisection(f , -5 , 0 , thresh, True)
		s_value , s_niter = rootfinder.Secant(f , -5 , 0 , thresh, True)
		nv , ni = rootfinder.Newton(f , fprime , -5 , thresh , True)
		th.append(np.log10(thresh))
		bi_iter.append(b_niter)
		n_iter.append(ni)
		s_iter.append(s_niter)
		thresh /=1.5
	print (nv , b_value , s_value)
	print (f(nv) , f(b_value) , f(s_value))
	plt.plot(th , bi_iter , label = "bisection")
	plt.plot(th , s_iter , label = "secant")
	plt.plot(th , n_iter , label = "Newton")
	plt.legend()
	plt.xlabel("log(threshhold)")
	plt.ylabel("Number of Iterations")
	plt.savefig("z.pdf")
	plt.close()
	
def problem_4():
	D = 1 * u.kpc
	a = 1 * u.au
	lamb = 21 * u.cm
	N_0 = .01 * u.pc / (u.cm ** 3)
	
	x = [-2 * u.au]
	xp = [lens(x[-1] , lamb , N_0 , D , a)]
	while x[-1].value < 2:
		x.append(x[-1]+ .01 *  u.au)
		xp.append(lens(x[-1] , lamb , N_0 , D , a))
	
	for i in range(len(x)):
		plt.plot([x[i].value , xp[i].value] , [1 , 0] , color = 'b' , linewidth = .5)
	plt.xlim(0 , 2)
	plt.savefig("4.pdf")
	plt.close()
	print (xp[i])
	
		
problem_3()
problem_4()
