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
		
data = np.random.poisson(3 , 1000)

T1 = MCMC.MCMC(P , data)

X = T1.M_H(Q , 100000 , [10])
N = range(len(X))


hx = []
xr = []
pr = []
for i in X:
	hx.append(i[0])



for i in range(max(hx)):
	xr.append(i)
	pr.append(poiss(3 , i) * len(hx))
	
plt.hist(hx , bins = max(hx))
plt.plot(xr , pr)
plt.show()
