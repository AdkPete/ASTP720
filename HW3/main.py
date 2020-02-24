
import ode
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
import astropy.units as u
import astropy.constants as const


def stiff_ode(L):
	def nfunc(t , y):
		return [-1 * L * (y[0] - np.cos(t))]
	return nfunc
	
def stiff_sltn(t , L):
	A = (-1 * L**2 / (1 + L ** 2)) * np.exp(-L * t) + (L / (1 + L ** 2)) * np.sin(t)
	A += (L ** 2 / (1 + L ** 2)) * np.cos(t)
	return A

def pend(y, t, b, c):
     theta, omega = y
     dydt = [omega, -b*omega - c*np.sin(theta)]
     return dydt
     
def my_pend(b , c):
	def npend(t , y):
		theta , omega = y
		dydt = [omega, -b*omega - c*np.sin(theta)]
		return dydt
	return npend

def WDPressure(rho):

	'''
	Equation of State for a White Dwarf
	Takes in a density rho
	Returns a Pressure
	'''
	
	###mu_e = 2
	C = (1 / 20.0) * (3 / np.pi) ** (2.0 / 3)
	C *= const.h ** 2 / (const.m_e * const.u ** (5.0 / 3))
	P =  C * ((rho / 2.0) ** (5.0 / 3))
	return P.to(u.Ba)
	
	
def problem_2():
	b = 0.25
	c = 5.0
	y0 = [np.pi - 0.1, 0.0]
	t = np.linspace(0, 10, 101)
	print (pend(y0 , t[0] , b , c))
	sol = odeint(pend, y0, t, args=(b, c))
	
	###Now we use my codes
	A = ode.solve_ode(my_pend(b , c) , 0 , y0 , (t[1] - t[0]))
	mt , my = A.Heun(t[-1])
	mt , my = A.Forward_Euler(t[-1])
	mt , my = A.RK4(t[-1])
	print (sol[0])
	pt = []
	ptheta = []
	for i in range(len(t)):
		pt.append(t[i])
		ptheta.append(sol[i][0])
		
	pmy =[]
	for i in my:
		pmy.append(i[0])
	plt.plot(pt , ptheta , label = "scipy")
	plt.plot(mt , pmy , label = "FE")
	plt.legend()
	plt.show()
	

def problem_3(L):

	
	#L = 1
	stiff_solve = ode.solve_ode(stiff_ode(L) , 0 , [0] , 1e-2)
	tend = 5
	t_h , y_h = stiff_solve.Heun(tend)
	t_fe , y_fe = stiff_solve.Forward_Euler(tend)
	t_rk4 , y_rk4 = stiff_solve.RK4(tend)
	rt = []
	ry = []
	for i in range(len(t_h)):
		rt.append(t_h[i])
		ry.append(stiff_sltn(t_h[i] , L))
		
	
	f , ax = plt.subplots( 2 , 2 , sharex = True)
	ax[0 , 0].plot(t_h , y_h , label = "numerical")
	ax[0 , 0].plot(rt , ry , label = "real")
	ax[0 , 0].legend()
	ax[0 , 0].set_title("Heun's Method")
	
	ax[1 , 0].plot(t_fe , y_fe , label = "numerical")
	ax[1 , 0].plot(rt , ry , label = "real")
	ax[1 , 0].legend()
	ax[1 , 0].set_title("Forward Euler Method")
	
	ax[0 , 1].plot(t_rk4 , y_rk4 , label = "numerical")
	ax[0 , 1].plot(rt , ry , label = "real")
	ax[0 , 1].legend()
	ax[0 , 1].set_title("RK4 Method")
	plt.show()
	
def problem_4():
	
	###white dwarfs
	
	rho_c = 1e4 * (u.g / (u.cm ** 3))
	
	## we will solve in terms of P and M
		
	P_c = WDPressure(rho_c).value() ###Central Pressure in Ba
	
	
#problem_2()
#problem_3(150)
problem_4()
