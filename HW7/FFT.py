import numpy as np
import matplotlib.pyplot as plt
import sys

def slow_dft(x):

	'''
	Takes in an array x
	returns the dft of x
	this is a slow, brute force approach
	'''
	N = len(x)
	
	X = []
	
	for k in range(len(x)):
	
		V = 0
		
		for n in range(len(x)):
			
			V +=  x[n] * np.exp(-2j *  np.pi * k * n / N)
		
		X.append(V)
		
	return np.array(X)
	

def slow_ift(x):

	'''
	Takes in an array x
	returns the inverse dft of x
	this is a slow, brute force approach
	'''
	N = len(x)
	
	X = []
	
	for k in range(len(x)):
	
		V = 0
		
		for n in range(len(x)):
			
			V +=  x[n] * np.exp(2j *  np.pi * k * n / N)
		
		X.append(V / N)
		
	return np.array(X)

	
def hpf(H , f , cut):

	'''
	This is a high pass filter in frequency space
	sets H to 0 at all frequencies below cut
	'''
	
	for i in range(len(H)):
		if f[i] > cut:
			H[i] = 0
	return H
	
def lpf(H , f , cut):

	'''
	This is a low pass filter in frequency space
	sets H to 0 at all frequencies above cut
	'''
	for i in range(len(H)):
		if f[i] < cut:
			H[i] = 0
	return H


	
def fft(x):

	'''
	This function will compute a fast fourer transform using the Cooley-Tukey Algorithm
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

def ifft(x):

	'''
	This function will compute an inverse fast fourer transform using the Cooley-Tukey Algorithm
	Takes in an array x
	returns the DFT of x
	'''



	N = len(x)

	w = np.exp(2j * np.pi / N)

	if N % 2 != 0:
		print ("Warning , N not a power of 2.")
		sys.exit()
		
	elif N <= 16:
		return slow_ift(x) / N
		
	else:
		
		##First we calculate the dfft at the even indices

		E = ifft(x[::2])
		
		##Now for the odd indices
		O = ifft(x[1::2])
		
		
		wk = []
		for k in range(N):
			wk.append( np.exp(2j * np.pi * k / N))
		
		wk = np.array(wk)


		X1 = E + wk[:int(N / 2)] * O
		X2 = E + wk[int(N / 2):] * O
		
		return list(X1) + list(X2)

def log_data(x , T):
	
	'''

	Takes in your original data set (time domain)
	T is the total time for which you were observing
	returns log10 of the fourier transform, and log10 of the frequency
	'''
	
	X = np.array(fft(x))
	X = np.log10(X)
	
	f = []
	for k in range(1 , len(X)):
		f.append(k / T)
	f = np.log10(np.array(f))
	
	X = np.delete(X , 0)
	return X , f
	
def filter(x , T , lf , hf):
	'''
	takes in a data set x, observed over a time T.
	filters out everything with f > hf and f < lf
	returns the inverse fourier transform of this data set
	'''
	
	
	X = np.array(np.fft.fft(x))
	
	f = []
	for k in range(0 , len(X)):
		f.append(k / T)
		
	X = hpf(X , f , hf)
	X = lpf(X , f , lf)

	return ifft(X)
	
	
	

	
