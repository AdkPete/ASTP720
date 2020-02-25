
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

###r is a radius in cm

	Pressure = y[0]
	Mass = y[1]
	if Pressure < 0:
		return [0 , 0]
	G = (const.G.to(u.cm ** 3 / (u.g * u.s ** 2))).value
	rho = NSrho((Pressure * u.Ba).to(u.GPa)).to(u.g / (u.cm ** 3)).value
	c = const.c.to(u.cm / (u.s)).value
	
	if Mass == 0:
		dp_dr = 0
	else:
		dp_dr = (-1 * G * Mass * rho / (r ** 2))
		dp_dr *= (1 + Pressure / (rho * c * c))
		dp_dr *= (1 + 4 * np.pi * r ** 3 * Pressure / (Mass * c * c))
		dp_dr *= (1 - 2 * G * Mass / (r * c * c)) ** -1
	dm_dr = 4 * np.pi * r ** 2 * rho

	return [dp_dr , dm_dr]
	
def problem_2():
	b = 0.25
	c = 5.0
	y0 = [np.pi - 0.1, 0.0]
	t = np.linspace(0, 10, 101)
	print (pend(y0 , t[0] , b , c))
	sol = odeint(pend, y0, t, args=(b, c))
	
	###Now we use my codes
	A = ode.solve_ode(my_pend(b , c) , 0 , y0 , (t[1] - t[0]))
	t_heun , y_heun = A.Heun(t[-1])
	t_fe , y_fe = A.Forward_Euler(t[-1])
	t_rk4 , y_rk4 = A.RK4(t[-1])
	print (sol[0])
	pt = []
	ptheta = []
	for i in range(len(t)):
		pt.append(t[i])
		ptheta.append(sol[i][0])
		
	fe = []
	for i in y_fe:
		fe.append(i[0])
		
	heun = []
	for i in y_heun:
		heun.append(i[0])
	
	rk4 = []
	for i in y_rk4:
		rk4.append(i[0])
	f , ax = plt.subplots(3 , 1 , sharex = True)
	ax[0].plot(pt , ptheta , label = "scipy")
	ax[0].plot(t_fe , fe , label = "Forward Euler")
	#ax[0].set_xlabel("time")
	ax[0].set_ylabel("theta")
	ax[0].legend()
	
	ax[1].plot(pt , ptheta , label = "scipy")
	ax[1].plot(t_heun , heun , label = "Heun")
	#ax[1].set_xlabel("time")
	ax[1].set_ylabel("theta")
	ax[1].legend()
	
	ax[2].plot(pt , ptheta , label = "scipy")
	ax[2].plot(t_rk4, rk4 , label = "RK4")
	ax[2].set_xlabel("time")
	ax[2].set_ylabel("theta")
	ax[2].legend()

	plt.show()
	
def stiff(L):
	
	#L = 1
	h = 1e-2
	stiff_solve = ode.solve_ode(stiff_ode(L) , 0 , [0] , h)
	tend = 5
	t_h , y_h = stiff_solve.Heun(tend)
	t_fe , y_fe = stiff_solve.Forward_Euler(tend)
	t_rk4 , y_rk4 = stiff_solve.RK4(tend)
	rt = []
	ry = []
	for i in range(len(t_h)):
		rt.append(t_h[i])
		ry.append(stiff_sltn(t_h[i] , L))
		
	heun = []
	rk4 = []
	fe = []
	
	t = 0
	
	for i in y_h:
		heun.append(abs(i - stiff_sltn(t , L)))
		t += h
	
	t = 0
	for i in y_fe:
		fe.append(abs(i - stiff_sltn(t , L)))
		t += h
	t = 0
	for i in y_rk4:
		rk4.append(abs(i - stiff_sltn(t , L)))
		t += h
		
	return max(fe) , max(heun) , max(rk4)
	
	'''
	f , ax = plt.subplots( 1 , 3 , sharex = True)
	ax[0].plot(t_h , heun , label = "numerical")
	#ax[0].plot(rt , ry , label = "real")
	#ax[0].legend()
	ax[0].set_title("Heun's Method")
	
	ax[1].plot(t_fe , fe , label = "numerical")
	#ax[1].plot(rt , ry , label = "real")
	#ax[1].legend()
	ax[1].set_title("Forward Euler Method")
	
	ax[2].plot(t_rk4 , rk4 , label = "numerical")
	#ax[2].plot(rt , ry , label = "real")
	#ax[2].legend()
	ax[2].set_title("RK4 Method")
	plt.show()
	'''
	
def problem_3():

	L = 1
	fe = []
	h = []
	Ls = []
	rk = []
	while L < 100:
		a , b  ,c = stiff(L)
		fe.append(a)
		h.append(b)
		rk.append(c)
		Ls.append(L)
		L += 1
	

	f , ax = plt.subplots( 3 , 1 , sharex = True)
	
	ax[0].plot(Ls , fe)
	ax[0].set_title("Forward Euler Method")
	
	ax[1].plot(Ls , h)
	ax[1].set_title("Heun's Method")
	
	ax[2].plot(Ls , rk)
	ax[2].set_title("RK4 Method")
	plt.show()
	
		
		
	
	
def problem_4():
	
	###white dwarfs
	
	rc = 1e4
	M_total = []
	R = []
	h = 5e6
	while rc <= 1e6:
	
	
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
	Pc = NSPressure(rc)
	Rad = []
	M_total = []
	while rc.value < 1e17:
		P_c = NSPressure(rc).value ###Central Pressure in Ba
		
		Neutron_Star = ode.solve_ode(NS_func , 1 , [P_c , 0] , h)
		
		rend = 5e9
		nsr , nsy = Neutron_Star.Forward_Euler(rend)

		for i in range(len(nsr)):
			if nsy[i][0] < 0:
				print (nsy[i])
				R = nsr[i]
				M = nsy[i][1]
				Rad.append((R * u.cm).to(u.km).value)
				M_total.append( (M * u.g).to(u.M_sun).value)
				break
		rc *= 1.2
	plt.plot(M_total , Rad)
	plt.xlabel("Mass (M_sun)")
	plt.ylabel("Radius (km)")
	plt.show()
problem_2()
problem_3()
problem_4()
problem_5()
