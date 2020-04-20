
import matplotlib.pyplot as plt
import astropy.units as u
import astropy.constants as const

import numpy as np

M = 908 * 1e12 * u.M_sun
R = 2 * u.Mpc

def F(X , Y):
	X *= u.Mpc
	Y *= u.Mpc
	
	zp = 0 * u.Mpc
	
	X -= 7 * u.Mpc
	Y -= 5 * u.Mpc
	
	def interior(x , y):
		return -1 * const.G * M * (3 * R ** 2 - (x ** 2 + y ** 2)) / (2 * R ** 3).to(u.Mpc ** 2 / u.Gyr ** 2)
	def exterior(x , y):
		return -1 * const.G * M / (np.sqrt(X ** 2 + Y ** 2))
	#Z = np.piecewise((X , Y), [np.sqrt(X ** 2 + Y ** 2) < R] , [interior , exterior])
	Z = np.select([ np.sqrt(X** 2 + Y ** 2) < R , np.sqrt(X** 2 + Y ** 2) >= R ] , [(-1 * const.G * M / (np.sqrt(X ** 2 + Y ** 2))).to(u.Mpc ** 2 / u.Gyr ** 2) , (-1 * const.G * M * (3 * R ** 2 - (X ** 2 + Y ** 2)) / (2 * R ** 3)).to(u.Mpc ** 2 / u.Gyr ** 2)])
	

	return Z
X = np.linspace(-30 , 30 , 500)
Y = np.linspace(-30 , 30 , 500)

X , Y = np.meshgrid(X , Y)


Z = F(X , Y)

plt.contour(X , Y , Z)
plt.xlim(-30 , 30)
plt.ylim(-30 , 30)
plt.show()
		
