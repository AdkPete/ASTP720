import nbody
import astropy.units as u
import time , sys
import numpy as np

F1 = np.load("galaxies0.npy")
F2 = np.load("galaxies1.npy")

if sys.argv[-1] == "1":
	Res = nbody.restart()
else:
	nbody.set_params("exp.txt")
		
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
		
	nx = 1 + np.random.rand() / 10.
	ny = (5 + np.random.rand() / 10. )
	nz = (5 + np.random.rand() / 10. )
	
	test_particle = nbody.Particle(nx * u.Mpc , ny * u.Mpc , nz * u.Mpc , 1e12 * u.M_sun , id)
	test_particle.r.append(nbody.	Vector([nx * u.Mpc , ny * u.Mpc , nz * u.Mpc]))
	test_id = id
	print (test_id)
	IC.append(test_particle)
	IC = nbody.Sim(IC)
	
	N = nbody.Barnes_Hut(IC)

	
