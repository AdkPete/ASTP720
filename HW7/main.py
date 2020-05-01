import numpy as np
import matplotlib.pyplot as plt
import FFT

def transform():
	'''
	Loads data from our data file
	runs a fourier transform
	saves the transformed data
	'''
	
	h = np.load("strain.npy")
	T = len(h) * 60 ###Total observation time in seconds

	Hk , freq = FFT.log_data(h , T)
	
	np.save("Transform" , [ Hk , freq ])
	
	mkplot()

def mkplot():

	'''
	Produces a log-log plot of our transformed data
	'''
	
	A = np.load("Transform.npy")
	
	freq = A[1]
	Hk = A[0]
	

	plt.plot(freq , Hk)
	plt.xlabel("log(Frequency (Hz))")
	plt.ylabel("log(Hk)")
	#plt.yscale("log")
	#plt.xscale("log")
	plt.xlim(-8 ,  -2)
	plt.savefig("Figure1.pdf")
	plt.close()
	
def measurements():

	'''
	This function reads in our data
	Prints out estimates for R and M
	Also prints the frequency of our signal
	And the amplitude of our signal
	'''
	
	D = 12
	
	A = np.load("Transform.npy")
	
	freq = A[1]
	Hk = A[0]
	
	rk = []
	rf = []
	
	for i in range(len(freq)):
		if freq[i] > -3 and freq[i] < -2:
			rf.append(freq[i])
			rk.append(Hk[i])
			
	peakh = max(rk)
	
	f = rf[rk.index(peakh)]
	

	f = 10 ** f.real
	
	print ("The frequency for or GW signal is {} Hz".format(f))

	amplitude = 10 ** peakh 
	
	
	amp = (np.sqrt(amplitude.real ** 2 + amplitude.imag ** 2) / (2 ** 19))
	
	print ("The amplitude for or GW signal is {}".format(amp))
	
	R = (amp / (2.6e-21)) * (1e-4 / f ) ** 4 * (D)
	R = R ** (1./5)
	
	print ("The seperation is R = {} Solar Radii".format(R))
	
	
	
	M = (R ** 1.5 * f / 1e-4) ** 2
	
	print ("The Mass is M = {} Solar Masses".format(M))
	
	
def run_filter():
	h = np.load("strain.npy")
	T = len(h) * 60 ###Total observation time in seconds
	
	X = FFT.filter(h , T , 1e-3 , 1e-2)
	t = []
	for i in range(len(X)):
		t.append(i / (60 * 24)) ##time in days
	
	plt.plot(t , np.array(X) / 1e-18)
	plt.xlabel("Time (days)")
	plt.ylabel("Strain / 1e-18 ")
	plt.savefig("Figure2.pdf")
	plt.close()
	
	
	plt.plot(t , np.array(h)/1e-18 , color = "green")
	plt.xlabel("Time (days)")
	plt.ylabel("Strain / 1e-18 ")
	plt.savefig("Figure3.pdf")
	
#transform()
measurements()
mkplot()
#run_filter()

