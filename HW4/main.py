import numpy as np
import nbody
import astropy.units as u

F1 = np.load("galaxies0.npy")
F2 = np.load("galaxies1.npy")

IC = []
for i in F1:
	IC.append(nbody.Particle(i[0] * u.Mpc , i[1] * u.Mpc , i[2] * u.Mpc , 1e12 * u.M_sun))
	
	