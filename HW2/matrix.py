import unittest
import numpy as np
import copy

class Matrix:
	
	def __init__(self , elements , size):
		self.size = size
		self.elements = elements
		
	def get(self , i , j):
		'''
		One indexed method to get a particular element
		The one indexing is useful since many matrix equations are one indexed.
		'''
		return self.elements[i - 1][j - 1]
		
	def __eq__(self , other):
		
		if self.size != other.size:
			return False
		for row in range(len(self.elements)):
			for col in range(len(self.elements[row])):
				if self.elements[row][col] != other.elements[row][col]:
					return False
		return True
		
	def print(self):
		for row in range(len(self.elements)):
			row_out = ""
			for col in range(len(self.elements[row])):
				row_out += str(self.elements[row][col]) + "\t"
			print (row_out)
		
	def __add__(self , other):
	
		New = copy.copy(Matrix(self.elements  , self.size))
		if self.size != other.size:
			 raise SystemExit('Attempted to add matrices that are not the same size')
		for row in range(len(self.elements)):
			for col in range(len(self.elements[row])):
				New.elements[row][col] += other.elements[row][col]
		return New

class TestMatrix(unittest.TestCase):
	def test_get(self):
		A = Matrix([[0 , 1 ] , [2 , 3]] , (2 , 2))
		self.assertEqual(A.get(1,1) , 0)
		
	def test_add(self):
		A = Matrix([[0 , 1 , 2 ] , [0 , 1 , 2 ] , [0 , 1 , 2 ]] , (3 , 3) )
		B = copy.copy(A)
		Answer = Matrix([[0 , 2 , 4] , [0 , 2 , 4] , [0 , 2 , 4]] , (3 , 3))
		C = A + B
		self.assertEqual(C , Answer)
		
		
if __name__ == "__main__":
	unittest.main()
