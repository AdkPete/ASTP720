from least_squares import *
import numpy as np
import plotting as pl
import scipy.optimize as opt
from normal import normal_method

def read_file(fname):

	'''
	reads in our data set
	takes in a file name called fname
	returns log(P) , M , Z, each of which being an array
	'''
	
	f = open(fname)
	Z = []
	M = []
	P = []
	for i in f.readlines():
		
		if i[0] == "#":
			continue
		else:
			P.append(np.log10(float(i.split(",")[1])))
			m = float(i.split(",")[3])
			distance = float(i.split(",")[2])
			d = distance = 1000.0 ##pc
			E = float(i.split(",")[-2])
			Rv = 3.1
			Av = Rv * E
			
			AM = m - 5 * np.log10(d) + 5 + Av
			M.append(AM)
			Z.append(float(i.split(",")[-1]))
	return P , M , Z

def fit_surface(alpha , beta ,  gamma):

	
	def f(P , Z):
		return alpha + beta * P + gamma * Z
	return f
		


def eval(P , M , Z):

	'''
	I used this for testing
	returns a function to evaluate a fit with the given data
	'''
	
	def f(x):
		alpha = x[0]
		beta = x[1]
		gamma = x[2]
		v = 0
		
		for i in range(len(P)):
			v += M[i] - (alpha + beta * P[i] + gamma * Z[i])
		
		return abs(v)
	
	return f
	



LP , M , Z = read_file("cepheid_data.txt")




alpha , beta , gamma , cov = more_general_fit(M , LP , Z)

#mk_3d(M , LP , Z , alpha , beta ,
f1 = fit_surface(alpha , beta , gamma)
pl.data_3d_with_surface(LP , Z , M , f1)
print (alpha , beta , gamma , "alpha")
print (alpha + beta * np.log10(11.36) + -0.080 * gamma)
f = eval(LP , M , Z)
print ("Our errors are sa = {} , sb = {} , sg = {}".format(np.sqrt(cov.elements[0][0]) , np.sqrt(cov.elements[1][1]) , np.sqrt(cov.elements[2][2])))
print (min(Z) , max(Z) , np.mean(Z))
cplot(LP , M , alpha , beta , gamma)


