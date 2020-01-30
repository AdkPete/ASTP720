import rootfinder
import interpolation
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
from astropy import constants as const

def	f(t):
	'''
	This function returns N_e(t) for the pseudo-isothermal sphere
	Takes in a number t which is equal to x / r_c
	'''
	return (1 + t ** 2) ** (-0.5) - .5

def fprime(t):
	'''
	This function returns the derivative of N_e(t) with respect to t
	This is meant to be used with the Newton method
	As above, t = x / r_c
	'''
	
	return -1 * t / ((1 + t ** 2) ** 1.5)
	
def lens(x , wlength , N0 , D , a):

	'''
	This is the lense equation given for the gaussian lens.
	Takes in a position, and several parameters of the lense.
	returns x prime
	'''
	
	#re = const.e.esu ** 2 / (const.m_e * const.c ** 2)
	re = 2.817940328e-15 * u.m
	C = ((wlength ** 2 * N0 * re * D) /(np.pi * a * a)) * np.exp(-((x / a) ** 2))
	return x * (1 + C)
	
def L2(x , wlength , N0 , D , a , rc):
	'''
	This is the second lense equation, which assumes that the lense is an isothermal sphere.
	Takes in a position, and several parameters of the lense.
	returns x prime
	'''
	re = 2.817940328e-15 * u.m
	C = ((wlength ** 2 * N0 * re * D) /(np.pi * a * a)) * ((1 + (x / rc) ** 2) ** (-0.5))
	return x * (1 + C)
	
def problem_3():
	
	'''
	My code for problem 3. This will determine the FWHM and produce the desired plots.
	'''
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
	plt.plot(th , bi_iter , label = "Bisection")
	plt.plot(th , s_iter , label = "Secant")
	plt.plot(th , n_iter , label = "Newton")
	plt.legend()
	plt.xlabel("log(threshold)")
	plt.ylabel("Number of Iterations")
	plt.show()
	
	print (nv , b_value , s_value)
	
def problem_4():
	D = 1 * u.kpc
	a = 1 * u.au
	lamb = 21 * u.cm
	N_0 = .01 * u.pc / (u.cm ** 3)
	rc = 1 * u.au
	
	x = [-10 * u.au]
	xp = [lens(x[-1] , lamb , N_0 , D , a)]
	xp2 = [L2(x[-1] , lamb , N_0 , D , a , rc)]
	while x[-1].value < 10:
		x.append(x[-1]+ .08 *  u.au)
		
		xp.append(lens(x[-1] , lamb , N_0 , D , a))
		xp2.append(L2(x[-1] , lamb , N_0 , D , a , rc))
	
	for i in range(len(x)):
		plt.plot([x[i].value , xp[i].value] , [1 , 0] , color = 'b' , linewidth = .5)
	#plt.xlim(-5 , 5)
	plt.xlabel("x (AU)")
	plt.ylabel("Distance (Kpc)")
	plt.show()
	
	for i in range(len(x)):
		plt.plot([x[i].value , xp2[i].value] , [1 , 0] , color = 'b' , linewidth = .5)
	#plt.xlim(0 , 5)
	plt.xlabel("x (AU)")
	plt.ylabel("Distance (Kpc)")
	plt.show()
	
	
		
problem_3()
problem_4()
