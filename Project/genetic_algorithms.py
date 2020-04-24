import numpy as np

def float_to_bin(N , p):
	N = round(N , p)
	
	i = int(str(N).split(".")[0])
	f = int(str(N).split(".")[1])
	
	return i , f
	
def bin_to_float(N):
	
	if float(N) == 0:
		return 0
	i = str(N).split(".")[0]
	
	
	a = len(str(i)) - 1
	
	v = 0
	k = 0
	i = str(i)
	
	while a >= 0:
		v += int(i[k]) * (2 ** a)
		a -= 1
		k += 1
		
	if len(str(N).split(".")) == 1:
		return v
	
	f  = str(N).split(".")[1]
	a = len(str(f)) - 1
	
	vf = 0
	k = 0
	f = str(f)
	
	while a >= 0:
		vf += int(f[k]) * (2 ** a)
		a -= 1
		k += 1
		
	res = float(str(v) + "." + str(vf))
	

	return res
	
class Genetic:

	def __init__(self , f , precision , creatures = None):
		'''
		f is a function that takes in an array of N parameters and returns some form of goodness of fit / likelihood or similar.
		N_param is the number of parameters to be fit.
		'''
		
		self.f = f
		self.prec = precision
		
		self.creatures = creatures
		
	
		
		
class creature:
	def __init__(self , dna):
		self.dna = dna
		self.fitness = None
		
		
	



print (v)
