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
	
	
def gausslens(wlength , N0 , D , a):
	'''
	This is the lense equation given for the gaussian lens.
	Takes in several parameters of the lense.
	wlength is a wavelength, N0 is the central column density, a is a characteristic size
	D is the distance to the lense
	returns a function of x that rturns x prime
	'''
	
	#re = const.e.esu ** 2 / (const.m_e * const.c ** 2)
	re = 2.817940328e-15 * u.m
	return lambda x: x * (1 + ((wlength ** 2 * N0 * re * D) /(np.pi * a * a)) * np.exp(-((x / a) ** 2)))
	
def isothermallens(wlength , N0 , D , rc):
	'''
	This is the lense equation for the isothermal lens.
	Takes in several parameters of the lense.
	wlength is a wavelength, N0 is the central column density, rc is a characteristic size
	D is the distance to the lense
	returns a function of x that rturns x prime
	'''
	re = 2.817940328e-15 * u.m
	c = (((wlength ** 2 * re) / (2 * np.pi))) ###theta_r = c * d/dx N_e(x)
	theta_r = lambda x: c * (N0 * x * -1) / (rc ** 2 * (1 + (x / rc) ** 2) ** 1.5)
	

	return lambda x: (x + D * theta_r(x))
	
def ray_tracing(L , spacing , xlow , xhigh):
	'''
	Produces a ray tracing plot, like on the first page of the HW.
	L should be a function of x that returns x prime
	Such a function caan be produces using gausslens or isothermallens
	spacing is the space between incoming rays
	xlow and xhigh are the lowest and highest x values for the incoming rays
	No returns
	'''
	
	spacing *= u.au
	x = xlow * u.au
	
	
	while x.value < xhigh:
		plt.plot([x.value , L(x).value] , [1 , 0] , color = 'b' , linewidth = .5)
		
		x += spacing
		
	plt.xlabel("x (AU)")
	plt.ylabel("Distance (Kpc)")
	plt.show()
		
	
	
def problem_3():
	
	'''
	My code for problem 3. This will determine the FWHM and produce the desired plots.
	No returns
	Will produce a plot
	'''
	thresh = 1
	th = []
	bi_iter = []
	s_iter = []
	n_iter = []
	while thresh > 1e-10: ###loop generates data for our plots
	
		b_value , b_niter = rootfinder.bisection(f , -5 , 0 , thresh, True) #bisection method
		
		s_value , s_niter = rootfinder.Secant(f , -5 , 0 , thresh, True) #Newton method
		
		nv , ni = rootfinder.Newton(f , fprime , -5 , thresh , True) #Secant method
		
		th.append(np.log10(thresh))
		bi_iter.append(b_niter)
		n_iter.append(ni)
		s_iter.append(s_niter)
		
		thresh /=1.5 #lowers threshold
		
	###Now we plot the resutls
	plt.plot(th , bi_iter , label = "Bisection")
	plt.plot(th , s_iter , label = "Secant")
	plt.plot(th , n_iter , label = "Newton")
	plt.legend()
	plt.xlabel("log(threshold)")
	plt.ylabel("Number of Iterations")
	plt.show()
	
	###Print out the value of the root found
	print (nv , b_value , s_value)
	
def problem_4():
	D = 1 * u.kpc
	a = 1 * u.au
	lamb = 21 * u.cm
	N_0 = .01 * u.pc / (u.cm ** 3)
	rc = 1 * u.au
	
	GL = gausslens(lamb , N_0 , D , a)
	ray_tracing(GL , .1 , -8 , 8)
	
	IL = isothermallens(lamb , N_0 , D , rc)
	ray_tracing(IL , .075 , -4 , 4)
	
	
		
problem_3()
problem_4()
