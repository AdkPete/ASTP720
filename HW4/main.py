import numpy as np
import nbody
import astropy.units as u
import time

F1 = np.load("galaxies0.npy")
F2 = np.load("galaxies1.npy")

IC = []
for i in F1:
	IC.append(nbody.Particle(i[0] * u.Mpc , i[1] * u.Mpc , i[2] * u.Mpc , 1e12 * u.M_sun))
	

for i in range(len(F2)):
	x = F2[i][0] * u.Mpc
	y = F2[i][1] * u.Mpc
	z = F2[i][2] * u.Mpc
	IC[i].r.append(nbody.Vector([x , y , z]))


s = time.time()
N = nbody.Barnes_Hut(IC , 1000 * u.yr , t_end = 100 * u.yr)
print ("Step time is {}").format(str(time.time() - s))

