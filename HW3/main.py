
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
	
def WDrho(P):	
	'''
	Equation of State for a White Dwarf
	Takes in a Pressure P
	Returns a density
	'''
	C = (1 / 20.0) * (3 / np.pi) ** (2.0 / 3)
	C *= const.h ** 2 / (const.m_e * const.u ** (5.0 / 3))
	rho = (P / C) ** (3.0 / 5)
	rho *= 2
	return rho
	
def NSPressure(rho):

	'''
	Equation of State for a White Dwarf
	Takes in a density rho
	Returns a Pressure
	'''

	###mu_e = 2
	C = (1 / 20.0) * (3 / np.pi) ** (2.0 / 3)
	C *= const.h ** 2 / (const.m_n ** (8.0 / 3))
	P =  C * ((rho) ** (5.0 / 3))
	return P.to(u.Ba)

def NSrho(P):

	'''
	Equation of State for a White Dwarf
	Takes in a density rho
	Returns a Pressure
	'''

	C = (1 / 20.0) * (3 / np.pi) ** (2.0 / 3)
	C *= const.h ** 2 / (const.m_n ** (8.0 / 3))
	rho = (P / C) ** (3.0 / 5)
	return rho
	


	
def WD_func(r , y):

	###r is a radius in cm
	
	Pressure = y[0]
	Mass = y[1]
	if Pressure < 0:
		return [0 , 0]
	G = (const.G.to(u.cm ** 3 / (u.g * u.s ** 2))).value
	rho = WDrho((Pressure * u.Ba).to(u.GPa)).to(u.g / (u.cm ** 3)).value
	if Mass == 0:
		dp_dr = 0
	else:
		dp_dr = (-1 * G * Mass * rho / (r ** 2))
	dm_dr = 4 * np.pi * r ** 2 * rho
	
	return [dp_dr , dm_dr]
	
def NS_func(r , y):
	
	Pressure = y[0] * u.Ba
	Mass = y[1] * u.g
	
	r *= u.cm
	if Pressure < 0:
		return [0 , 0]
	G = (const.G.to(u.cm ** 3 / (u.g * u.s ** 2)))
	rho = NSrho((Pressure)).to(u.g / (u.cm ** 3))
	c = const.c.to(u.cm / (u.s))
	
	if Mass.value == 0:
		dm_dr = 4 * np.pi * r ** 2 * rho
		return [ 0 , dm_dr.value ]
	else:
		dp_dr = -1 * G * Mass * rho / (r ** 2)
		dp_dr *= (1 + Pressure / (rho * c ** 2))
		dp_dr *= (1 + 4 * np.pi * (r ** 3) * Pressure / (Mass * c ** 2))
		dp_dr *= (1 - 2 * G * Mass / (r * c ** 2)) ** (-1)
	
	dm_dr = 4 * np.pi * r ** 2 * rho
	
	return [ dp_dr.value , dm_dr.value]
	
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
	
	rc = 1e3
	M_total = []
	R = []
	h = 5e6
	while rc <= 1e7:
	
	
		rho_c = rc * (u.g / (u.cm ** 3))
		
		## we will solve in terms of P and M
			
		P_c = WDPressure(rho_c).value ###Central Pressure in Ba
		
		White_Dwarf = ode.solve_ode(WD_func , 1 , [P_c , 0] , h)
		
		rearth = 6.378e+8 ##cm
		r_end = 5 * rearth
		
		wdr , wdy = White_Dwarf.RK4(r_end)
		

		for i in range(len(wdr)):
			if wdy[i][0] < 0:
				R.append((wdr[i] * u.cm).to(u.R_sun).value)
				M_total.append(((wdy[i][1] * u.g).to(u.M_sun).value))
				break
		rc *= 1.5
	plt.plot(M_total , R)
	plt.xlabel("Mass (M_sun)")
	plt.ylabel("Radius (Solar Radii)")
	plt.show()
	
def problem_5():

	###neutron stars
	h = 1e4
	
	rc = 1e14 * u.g / (u.cm ** 3)
	Rad = []
	M_total = []
	while rc.value < 1e17:
		P_c = NSPressure(rc).value ###Central Pressure in Ba
		
		Neutron_Star = ode.solve_ode(NS_func , 1 , [P_c , 0] , h)
		
		rend = 5e6
		nsr , nsy = Neutron_Star.RK4(rend)

		for i in range(len(nsr)):
			if nsy[i][0] < 0:
				R = nsr[i]
				M = nsy[i][1]
				Rad.append((R * u.cm).to(u.km).value)
				M_total.append( (M * u.g).to(u.M_sun).value)
				break
		rc *= 2
	plt.plot(M_total , Rad)
	plt.xlabel("Mass (M_sun)")
	plt.ylabel("Radius (Solar Radii)")
	plt.show()
#problem_2()
#problem_3(150)
#problem_4()
problem_5()
