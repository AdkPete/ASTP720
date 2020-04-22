import MCMC
import numpy as np
import matplotlib.pyplot as plt

def poiss(L , k):
	return (L ** k * np.exp(-L)) / np.math.factorial(k)
def P(X , D):
	L = np.mean(D)
	k = X[0]
	
	return (L ** k * np.exp(-L)) / np.math.factorial(k)
	
def Q(X):
	a = np.random.rand()

	if X[0] == 0:
		return [X[0] + 1]
	if a > 0.5:
		return [X[0] + 1]
	else:
		return [X[0] - 1]
		
def mp(data):

	for i in range(max(hx)):
		xr.append(i)
		pr.append(poiss(lam , i) * len(hx))
	y,binEdges = np.histogram(data,bins=max(data))
	bincenters = range(max(data))
	width      = 0.25
	plt.bar(bincenters, y, width=width, color='b')
	plt.plot(xr , pr , color = 'r')
	plt.show()
		
lam = 4

data = np.random.poisson(lam , 1000)

T1 = MCMC.MCMC(P , data)

X = T1.M_H(Q , 1000000 , [10])[100::]
N = range(len(X))


hx = []
xr = []
pr = []
for i in X:
	hx.append(i[0])




print (np.mean(hx))

mp(hx)
