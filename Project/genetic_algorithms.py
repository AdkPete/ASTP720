import numpy as np
import sys

def bin_inc(string , ind):


	if ind > len(string):
		print ("error")
		sys.exit()
		
	if ind == 0:
	
		return "1" + string[1::]
	
	elif ind == len(string) - 1:
		return string[0:-1] + "1"
	
	else:
		return string[0:ind] + "1" + string[ind + 1::]

def float_to_bin(N , p):
	
	i = int(str(N).split(".")[0])
	
	dna_i = '0' * p[0]
	
	for k in range(len(dna_i)):

		if (2 ** (k)) > i:
			break
		elif i % (2 ** (k + 1)) != 0:
			dna_i = bin_inc(dna_i , len(dna_i) - 1 - k)
			i -= 2 ** k
			
	if len(str(N).split(".")) == 0:
		dna_f = '0' * p[1]
		return dna_i + "." + dna_f
	
	f = int(str(N).split(".")[1])
	
	dna_f = "0" * p[1]
	
	
	for k in range(len(dna_f)):

		if (2 ** (k)) > f:
			break
		elif f % (2 ** (k + 1)) != 0:
			dna_f = bin_inc(dna_f , len(dna_f) - 1 - k)
			f -= 2 ** k
	
	
	
	return str(dna_i) + "." + dna_f
	
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
		precision determines how large the dna strings will be
		should be array-like containint 2 elements.
		The first describes the number of values allowed before the decimal point.
		The second describes the number of values allowed after the decimal point.
		'''
		
		self.f = f
		self.prec = precision
		
		self.creatures = creatures
		
	
		
		
class creature:
	def __init__(self , dna):
		self.dna = dna
		self.fitness = None
		
		
	
N = 101.253	
bin = float_to_bin(15.74 , [8 , 8 ])
numb = bin_to_float(bin)
print (N)
