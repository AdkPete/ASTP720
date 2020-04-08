import copy
import numpy as np
class matrix:
	
	def __init__(self , elements):
		self.elements = np.array(elements)
			
	
	def __getitem__(self , i , j):
		return self.elements[i][j]
		
	def __setitem__(self , i , j , val):
		self.elements[i][j] = val
		
	def __mul__(self , other):
		return matrix(np.matmul(self.elements , other.elements))
		
	def transpose(self):
		N = matrix(self.elements.T)
		return N
		
	def inverse(self):
		N = matrix(np.linalg.inv(self.elements))
		return N

	def print(self):
		for i in self.elements:
			print (i)
