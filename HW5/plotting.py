import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

def data_3d_with_surface(x , y , z , f):
	'''
	x , y , and z should be arrays containing all of your data.
	f should be a function of x and y that produces a z value for the surface
	'''
	
	lx = np.linspace(min(x) , max(x))
	ly = np.linspace(min(y) , max(y))
	
	sx , sy = np.meshgrid(lx , ly)
	sz = f(sx , sy)
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x, y, z)
	ax.set_xlabel("Log(P)")
	ax.set_ylabel("[Fe/H]")
	ax.set_zlabel("M")
	ax.plot_wireframe(sx , sy , sz)
	plt.show()
	

