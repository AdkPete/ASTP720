import numpy as np
import matplotlib.pyplot as plt
import sys

def slow_dft(x):

	N = len(x)
	
	X = []
	
	for k in range(len(x)):
	
		V = 0
		
		for n in range(len(x)):
			
			V +=  x[n] * np.exp(-2j *  np.pi * k * n / N)
		
		X.append(V)
		
	return np.array(X)
	
def fft(x):

	'''
	This function will compute a fast fourir transform using the Cooley-Tukey Algorithm
	Takes in an array x
	returns the DFT of x
	'''
	
	
	
	N = len(x)
	
	w = np.exp(-2j * np.pi / N)
	
	if N % 2 != 0:
		print ("Warning , N not a power of 2.")
		sys.exit()
		
	elif N <= 16:
		return slow_dft(x)
		
	else:
		
		##First we calculate the dfft at the even indices

		E = fft(x[::2])
		
		##Now for the odd indices
		O = fft(x[1::2])
		
		
		wk = []
		for k in range(N):
			wk.append( np.exp(-2j * np.pi * k / N))
		
		wk = np.array(wk)


		X1 = E + wk[:int(N / 2)] * O
		X2 = E + wk[int(N / 2):] * O
		
		return list(X1) + list(X2)


