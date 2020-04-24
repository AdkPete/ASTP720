import numpy as np


class MCMC:
	
	def __init__(self , P , data):
	
		def f(X):
			#print (X)
			return P(X , data)
		self.P = P	
		self.f = f ##This is the desired posterior distribution
		self.data = data
		self.MAP = None
		self.MAPL = None
		
	
	def M_H(self , Q , N , X_0):
		'''
		This method will implement a Metropolis-Hastings algorithm.
		The posterior distribution will be self.P
		The proposal distribution is Q
		The number of steps to take is N
		X_0 is an array containing you starting location
		'''
	
		step = 0
		
		X = np.array([X_0])
		self.MAP = [ X_0 ]
		self.MAPL = self.f(X[0])
		
		while step < N:
			#print (X)
			
			Y = Q(X[-1])
			
			NL = self.f(Y)
			r =  NL / self.f(X[-1])
			if NL > self.MAPL:
				self.MAPL = NL
				self.MAP = Y
			if r >= 1:
				
				X = np.append(X , [ Y ] , axis = 0)
			else:
				U = np.random.rand()
				if U <= r:
					X = np.append(X , [ Y ] , axis = 0)
				
				else:
				
					X = np.append(X , [ X[-1] ] , axis = 0)
				
			
			step += 1
		return X , self.MAP , self.MAPL
