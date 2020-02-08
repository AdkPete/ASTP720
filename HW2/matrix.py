import unittest
import numpy as np

class Matrix:
	
	def __init__(self , size , elements = []):
		if elements == []:
			row = 0
			self.elements = []
			while row < size[0]:
				col = 0
				self.elements.append([])
				while col < size[1]:
					self.elements[-1].append(0)
					col += 1
				row += 1
		else:
			self.elements = elements
		self.size = size
		
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
		New = Matrix(self.size)
		if self.size != other.size:
			 raise SystemExit('Attempted to add matrices that are not the same size')
		for row in range(len(self.elements)):
			for col in range(len(self.elements[row])):
				New.elements[row][col] = self.elements[row][col] + other.elements[row][col]
		return New

	def transpose(self):
		New = Matrix((self.size[0] , self.size[1]))
		for row in range(len(self.elements)):
			for col in range(len(self.elements[row])):
				New.elements[col][row] = self.elements[row][col]
				
		return New
		
class TestMatrix(unittest.TestCase):
	
	
	def test_get(self):
		A = Matrix((2 , 2) , [[0 , 1 ] , [2 , 3]])
		self.assertEqual(A.get(1,1) , 0)
		
	def test_add(self):
		A = Matrix((3 , 3) , [[0 , 1 , 2 ] , [0 , 1 , 2 ] , [0 , 1 , 2 ]] )
		B = Matrix((3 , 3) , [[0 , 1 , 2 ] , [0 , 1 , 2 ] , [0 , 1 , 2 ]] )
		Answer = Matrix((3 , 3) , [[0 , 2 , 4] , [0 , 2 , 4] , [0 , 2 , 4]])
		C = A + B
		self.assertEqual(C , Answer)
		
	def test_transpose(self):
		A = Matrix( (3 , 3) , [[0 , 0 , 0 ] , [1 , 1 , 1 ] , [2 , 2 , 2 ]])
		AT = Matrix( (3 , 3) , [[0 , 1 , 2] , [0 , 1 , 2] , [0 , 1 , 2]])
		self.assertEqual(AT , A.transpose())
		
		
if __name__ == "__main__":
	unittest.main()
