import numpy as np
import nbody
import astropy.units as u
import time , sys

F1 = np.load("galaxies0.npy")
F2 = np.load("galaxies1.npy")

if sys.argv[-1] == "1":
	Res = nbody.restart()
else:

	if len(sys.argv) == 1:
		nbody.set_params("params.txt")
	elif len(sys.argv) >= 2:
		nbody.set_params(sys.argv[1])
		
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
	
	'''
	av_vel = nbody.avg_velocity(IC , 1000 * u.yr)
	print (av_vel)
	targ_distance = 2 * u.Mpc
	run_time = (targ_distance) / av_vel
	print ("Desired run time for a typical particle to move {} MPC is {} Myr.".format(targ_distance.value , run_time.to(u.Myr).value))
	'''


	N = nbody.change_dt(IC , 1000 * u.yr)
	s = time.time()
	#N = nbody.Barnes_Hut(IC)
	print ("Run time is {}".format(time.time() - s))
	
