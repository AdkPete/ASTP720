import numpy as np
import matplotlib.pyplot as plt
import FFT

def transform():
	h = np.load("strain.npy")
	T = len(h) * 60 ###Total observation time in seconds

	Hk , freq = FFT.log_data(h , T)

	np.save("Transform" , [ Hk , freq ])
	
	mkplot()

def mkplot():
	A = np.load("Transform.npy")
	
	freq = A[1]
	Hk = A[0]
	
	plt.plot(freq , Hk)
	plt.xlabel("log ( Frequency (Hz) )")
	plt.ylabel("log ( Hk ) ")
	plt.xlim(-8 , -2)
	plt.savefig("Figure1.pdf")
	
def measurements():
	
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
	print (f.real)
	f = 10 ** f.real
	amplitude = 10 ** peakh 
	amp = (np.sqrt(amplitude.real ** 2 + amplitude.imag ** 2) / (2 ** 19))
	
	R = (amp / (2.6e-21)) * (1e-4 / f ) ** 4 * (D)
	R = R ** (1/5)
	print (R)
	
	M = (R ** 1.5 * f / 1e-4) ** 2
	print (M)
	
	print (f , 1e-4 * np.sqrt(M) * R ** (-1.5))
	print (amp , 2.6e-21 * M * M * (1. / D) / R)
#transform()
measurements()
#mkplot()
