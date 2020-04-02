import numpy as np
import nbody
import astropy.units as u
import time

F1 = np.load("galaxies0.npy")
F2 = np.load("galaxies1.npy")
nbody.set_params("params.txt")
IC = []
id = 0
for i in F1:
	IC.append(nbody.Particle(i[0] * u.Mpc , i[1] * u.Mpc , i[2] * u.Mpc , 1e12 * u.M_sun , id))
	
	id += 1

for i in range(len(F2)):
	x = F2[i][0] * u.Mpc
	y = F2[i][1] * u.Mpc
	z = F2[i][2] * u.Mpc
	IC[i].r.append(nbody.Vector([x , y , z]))
	
IC = nbody.Sim(IC)


#N = nbody.change_dt(IC , 1000 * u.yr , 1e6 * u.yr , 3e6 * u.yr)
s = time.time()
N = nbody.Barnes_Hut(IC)
print ("Step time is {} seconds , {} minutes".format(str(time.time() - s) , str((time.time() - s) / 60.0)))

