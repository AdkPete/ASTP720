import numpy as np
import matplotlib.pyplot as plt

def slow_dft(x):

	N = len(x)
	
	X = []
	
	for k in range(len(x)):
	
		V = 0
		
		for n in range(len(x)):
			
			V +=  x[n] * np.exp(-2j *  np.pi * k * n / N)
		
		X.append(V)
		
	return X
	
def test_fft():
	x = np.linspace(0 , 10 , 100)
	t = np.linspace(0 , 10 , 100)
	x = np.sin(2 * np.pi * x)
	

	X = slow_dft(x)
	
	NX = np.fft.fft(x)
	
	plt.plot(t , X)
	plt.show()
	print (np.allclose(NX , X))


test_fft()
