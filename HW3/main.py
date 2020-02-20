
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

def problem_3(L):

	
	#L = 1
	stiff_solve = ode.solve_ode(stiff_ode(L) , 0 , [0] , 1e-5)
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
	
problem_3(10)
