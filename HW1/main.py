import rootfinder
import interpolation
import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u
def	f(t):
	return (1 + t ** 2) ** (-0.5) - .5

def fprime(t):
	return -1 * t / ((1 + t ** 2) ** 1.5)
	
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
	
def problem_4():
	D = 1 * u.kpc
	a = 1 * u.au
	lamb = 21 * u.cm
	N_0 = 0.01 * u.pc / (u.cm ** 3)
	
problem_3()
problem_4()
