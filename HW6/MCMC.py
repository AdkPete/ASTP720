import numpy as np


class MCMC:
	
	def __init__(self , P , data):
	
		def f(X):
			#print (X)
			return P(X , data)
		self.P = P	
		self.f = f ##This is the desired posterior distribution
		self.data = data
		
	
	def M_H(self , Q , N , X_0):
		'''
		This method will implement a Metropolis-Hastings algorithm.
		The posterior distribution will be self.P
		The proposal distribution is Q
		The number of steps to take is N
		X_0 is an array containing you starting location
		'''
	
		step = 0
		
		X = [X_0]
		
		while step < N:
			#print (X)
			
			Y = Q(X[-1])
			
			r = self.f(Y) / self.f(X[-1])
			
			if r >= 1:
				X.append(Y)
			else:
				U = np.random.rand()
				if U <= r:
					X.append(Y)
				
			
			step += 1
			
		return X
