
import lightkurve as lk
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt
import sys
import scipy.stats as stats
from genetic_algorithms import Genetic


def get_data():

	'''
	This function will read in a lightcurve
	Then it removes outliers and fills in any gaps
	Returns the frequencies and powers for a periodogram of the lightcurve
	'''
	
	search = lk.search_lightcurvefile('KIC10963065', cadence='short', mission='Kepler')
	#print(search)

	files = search[1:10].download_all()

	lc = files.PDCSAP_FLUX.stitch()

	lc = lc.remove_nans().remove_outliers().fill_gaps()

	lc = lc.to_periodogram(freq_unit=u.microHertz, maximum_frequency=3000, minimum_frequency=1500)

	frequency = lc.frequency
	
	power = lc.power
	

	return frequency , power
	

	
def peaks(frequency , power):

	'''
	
	This function will filter out most of the periodogram
	selects only the peaks caused by astrosiesmic oscillations.
	returns the filtered set of frequencies and powers
	'''
	
	
	nf = []
	npr = []
	
	cut = .00000075
	
	
	for i in range(len(power)):

		if power[i].value < cut or frequency[i] in nf:
			continue
		npr.append(power[i].value)
		nf.append(frequency[i].value)
	

	
	#return nf , npr
	fit_power = []
	fit_frequency = []
	
	npr = np.array(npr)
	nf = np.array(nf)
	
	while len(npr) > 0:
		
		max_pow = max(npr)
		
		index = np.argmax(npr)
		peak_freq = nf[index]
		
		fit_power.append(max_pow)
		fit_frequency.append(peak_freq)
		
		del_ind = []
		for i in range(len(nf)):
			if abs(nf[i] - peak_freq) < 10:
				del_ind.append(i)
				
		npr = np.delete(npr , del_ind)
		nf = np.delete(nf , del_ind)
		
		
	return fit_frequency , fit_power

	
def model(n , l , X):

	epsnl = 0
	
	D_Nu = X[0]
	d_nu = X[1]
	alpha = X[2]
	
	return D_Nu * (n + l / 2. + alpha + 1/4.) + epsnl
		
def least_squares(frequency , power):
	def f(X):
			
		mean = X[0]
		deviation = X[1]
		amp = X[2]
		
		lsq = 0
		for i in range(len(power)):
			
			lsq -= ( power[i] - amp * stats.norm.pdf(frequency[i] , mean , deviation) ) ** 2
			
		return lsq
		
	return f
	
def plot_fit(X , frequency , power):
	
	mean = X[0]
	deviation = X[1]
	amp = X[2] / 1000
	
	tf = []
	tfit = []
	for i in np.linspace(min(frequency) , max(frequency) , 1000):
	
		theory = amp * stats.norm.pdf(i , mean , deviation)
		tfit.append(theory)
		tf.append(i)
		
	plt.scatter(frequency , power)
	plt.plot(tf , tfit)
	plt.ylim(.00000008 , .0000017)
	plt.show()
	

		
#frequency , power = get_data()
#frequency , power = peaks(frequency , power)

#np.save("Temp" , [frequency , power])
'''
A = np.load("Temp.npy")
frequency = A[0]
power = A[1]
fit = least_squares(frequency , power)

X = [np.mean(frequency) , np.std(frequency) , 1.5]
print (fit(X))

#opt = Genetic(fit , [20 , 20] , [ [1000 , 4000] , [.01 , 1000 ] , [0 , .1 ] ] , 10000)
opt = Genetic(fit , [20 , 20] , [ [2100 , 2200] , [400 , 420 ] , [1 , 2 ] ] , 1000)
opt.update(10)

print (opt.creatures[0].get_params() , opt.creatures[0].fitness , fit(X) < opt.creatures[0].fitness , fit(X))
#print (X)
plot_fit(opt.creatures[0].get_params() , frequency , power)
'''

Teff = 6046
Nu_Max = 2197
dnu = 103

R = (135 / dnu) ** 2 * (Nu_Max / 3050) * (Teff / 5777) ** (1 / 2.)
M = (135 / dnu) ** 4 * (Nu_Max / 3050) ** 3 * (Teff / 5777) ** (3. / 2.)

print (R , M)

'''
plt.scatter(frequency , power)
plt.ylim(.0000008 , .0000017)
plt.show()
'''


