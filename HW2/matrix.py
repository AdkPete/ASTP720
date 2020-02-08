import unittest
import numpy as np

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

class TestMatrix(unittest.TestCase):
	def test_get(self):
		A = Matrix([[0 , 1 ] , [2 , 3]] , (2 , 2))
		self.assertEqual(A.get(1,1) , 0)
		
if __name__ == "__main__":
	unittest.main()
