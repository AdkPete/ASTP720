
import matplotlib.pyplot as plt
from matrix import matrix
from mpl_toolkits.mplot3d import Axes3D


def lin_fit(y , t , eps = None):
	###We are going to solve a matrix problem to find m and b, where y = m * t + b
	n = len(y)
	Y = matrix(y)
	if eps != None:
		Err = matrix(eps)
	X = []
	
	for i in range(n):
		X.append([1 , t[i]])
	
	X = matrix(X)
	
	Xdagger = X.transpose()
	XdagX = Xdagger * X
	
	theta = XdagX.inverse() * Xdagger
	theta *= Y

	return theta.elements[0] , theta.elements[1]
	
def more_general_fit(y , t1 , t2 , eps = None):
	###We are going to solve a matrix problem to find m and b, where y = m * t + b
	n = len(y)
	Y = matrix(y)
	X = []
	if eps != None:
		Err = matrix(eps)
	
	for i in range(n):
		X.append([1 , t1[i] , t2[i]])
	
	X = matrix(X)
	
	Xdagger = X.transpose()
	XdagX = Xdagger * X
	
	theta = XdagX.inverse() * Xdagger
	theta *= Y

	return theta.elements[0] , theta.elements[1] , theta.elements[2] , (XdagX.inverse())

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
			

def cplot(LP , M , alpha , beta , gamma):
	'''
	makes a plot of M vs P at some fixed Z value
	
	'''
	zl = [-1.5 , -1 , -0.5 , 0 , 0.5 , 1 , 1.5]
	for tz in zl:
	
		pfit = []
		mfit = []
		
		p = min(LP)
		while p < max(LP):
			pfit.append(p)
			mfit.append(alpha + beta * p + gamma * tz)
			p += 0.01
		plt.plot(pfit , mfit , label = "[Fe/H] = {}".format(tz))
		
	plt.scatter(LP , M)
	plt.xlabel("log(P)")
	plt.ylabel("M")
	plt.legend()
	plt.show()

		
