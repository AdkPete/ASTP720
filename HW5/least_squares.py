
import matplotlib.pyplot as plt
from matrix import matrix
from mpl_toolkits.mplot3d import Axes3D


def lin_fit(y , t):
	###We are going to solve a matrix problem to find m and b, where y = m * t + b
	n = len(y)
	Y = matrix(y)
	X = []
	
	for i in range(n):
		X.append([1 , t[i]])
	
	X = matrix(X)
	
	Xdagger = X.transpose()
	XdagX = Xdagger * X
	
	theta = XdagX.inverse() * Xdagger
	theta *= Y

	return theta.elements[0] , theta.elements[1]
	
def more_general_fit(y , t1 , t2):
	###We are going to solve a matrix problem to find m and b, where y = m * t + b
	n = len(y)
	Y = matrix(y)
	X = []
	
	for i in range(n):
		X.append([1 , t1[i] , t2[i]])
	
	X = matrix(X)
	
	Xdagger = X.transpose()
	XdagX = Xdagger * X
	
	theta = XdagX.inverse() * Xdagger
	theta *= Y

	return theta.elements[0] , theta.elements[1] , theta.elements[2]

def mk_lin_plot(y , t , m , b):
	
	start_t = min(t)
	end_t = max(t)
	h = (end_t - start_t) / 1000.0
	tf = start_t
	
	fit_t = []
	fit_y = []
	
	while tf < end_t:
		
		fit_t.append(tf)
		fit_y.append(m * tf + b)
		tf += h
		
	plt.plot(fit_t , fit_y)
	plt.scatter(t , y)
	plt.show()
			
def mgplot(y , t1 , t2 , a , b , c):
	
	start_t = min(t1)
	end_t = max(t1)
	h = (end_t - start_t) / 1000.0
	tf = start_t
	
	fit_t = []
	fit_y = []
	t2 = 1
	while tf < end_t:
		
		fit_t.append(tf)
		fit_y.append(a + b * tf + c * t2)
		tf += h
		
	plt.plot(fit_t , fit_y)
	plt.scatter(t1 , y)
	plt.show()

def mk_3d(z , x , y , a , b , c):
	
	
	fit_1x = []
	fit_1y = []
	fit_1z = []
	
	for i in range(0 * 100 , 2 * 100):
		nx = i / 100.0
		ny = 0
		fit_1y.append(ny)
		fit_1x.append(nx)
		fit_1z.append(a + b * nx + c * ny)
		
	fit_2x = []
	fit_2y = []
	fit_2z = []
	
	for i in range(0 * 100 , 2 * 100):
		nx = i / 100.0
		ny = 0.5
		fit_2y.append(ny)
		fit_2x.append(nx)
		fit_2z.append(a + b * nx + c * ny)
		
		
	'''
	plt.plot(fit_t , fit_y)
	plt.scatter(t1 , y)
	plt.show()
	'''
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.scatter(x , y , z , color = "g")
	ax.plot(fit_1x , fit_1y , fit_1z)
	ax.plot(fit_2x , fit_2y , fit_2z)
	plt.show()
