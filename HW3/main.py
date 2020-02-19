
import ode
import matplotlib.pyplot as plt
import numpy as np


def stiff_ode(L):
	def nfunc(t , y):
		return [-1 * L * (y[0] - np.cos(t))]
	return nfunc
	
def stiff_sltn(t , L):
	A = (-1 * L**2 / (1 + L ** 2)) * np.exp(-L * t) + (L / (1 + L ** 2)) * np.sin(t)
	A += (L ** 2 / (1 + L ** 2)) * np.cos(t)
	return A

def problem_3():

	
	L = 100
	stiff_solve = ode.solve_ode(stiff_ode(L) , 0 , [0] , 1e-5)
	tend = 5
	t_100 , y_100 = stiff_solve.Heun(tend)
	
	rt_100 = []
	ry_100 = []
	for i in range(len(t_100)):
		rt_100.append(t_100[i])
		ry_100.append(stiff_sltn(t_100[i] , L))
		
	plt.plot(t_100 , y_100 , label = "numerical")
	plt.plot(rt_100 , ry_100 , label = "real")
	plt.legend()
	plt.show()
	
problem_3()
