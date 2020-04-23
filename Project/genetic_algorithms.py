import numpy as np

class Genetic:

	def __init__(self , f , param_string , creatures = None):
		'''
		f is a function that takes in an array of N parameters and returns some form of goodness of fit / likelihood or similar.
		N_param is the number of parameters to be fit.
		'''
		
		self.f = f
		self.PS = param_string
		
		self.creatures = creatures
		
	def dna_to_params(self , dna):
		
		X = np.array( [] )
		pos = 0
		for i in self.PS.split("--"):
			A = 0
			A += 8 * int(dna[pos]) + 4 * int(dna[pos + 1]) + 2 * int(dna[pos + 2]) + int(dna[pos + 3])
			pos += 4
			
			dp = 1
			while dp <= int(i.split(".")[1]):
				A += 10 ** (-dp) * (8 * int(dna[pos]) + 4 * int(dna[pos + 1]) + 2 * int(dna[pos + 2]) + int(dna[pos + 3]))
				pos += 4
				dp += 1
			
			
			X = np.append(X , A * (10 ** (int(i.split(".")[-1]))))
		
		return X
		
class creature:
	def __init__(self , dna):
		self.dna = dna

def test():
	A = Genetic(1 , "1.2.0")
	#A.dna_to_params(5)
	B =  "101101001100"
	
	print (A.dna_to_params(B))
	
test()
