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
		
	def __mul__(self , other):
		N = Matrix((self.size[0] , other.size[1]))
		
		for row in range(len(N.elements)):
			for col in range(len(N.elements[row])):
				V = 0
				for i in range(len(self.elements)):
					V += self.elements[row][i] * other.elements[i][col]
				
				N.elements[row][col] = V
		
		return N
		
	def trace(self):
		TR = 0
		if self.size[0] != self.size[1]:
			raise SystemExit('Cannot calculate the trace of a non-square matrix')
		for i in range(self.size[0]):
			TR += self.elements[i][i]
		
		return TR
		
	def residual(self , i , j):
		Res = Matrix((self.size[0] - 1 , self.size[1] - 1))
		Res.elements = []
		for row in range(len(self.elements)):
			if row == i:
				continue
			Res.elements.append([])
			for col in range(len(self.elements)):
				if col == j:
					continue
				Res.elements[-1].append(self.elements[row][col])
				
		return Res
					
				
	def determinant(self):
		if self.size == (2 , 2):
			return self.elements[0][0] * self.elements[1][1] - self.elements[0][1] * self.elements[1][0]
		###Now we handle larger matrices
		det = 0
		for i in range(len(self.elements)):
			Res = self.residual(0 , i)
			det += (-1) ** (i + 2) * self.elements[0][i] * Res.determinant()
		
		return det
		
	def inverse(self):
		New = Matrix(self.size)
		D = self.determinant()
		for row in range(len(self.elements)):
			for col in range(len(self.elements[row])):
				Cij = (-1) ** (row + col + 2) * self.residual(row , col).determinant()
				New.elements[row][col] = Cij / D
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
		
	def test_mult(self):
		A =  Matrix((3 , 3) , [[0 , 1 , 2 ] , [0 , 1 , 2 ] , [0 , 1 , 2 ]] )
		B =  Matrix( (3 , 3) , [[0 , 0 , 0 ] , [1 , 1 , 1 ] , [2 , 2 , 2 ]])
		AB = Matrix( (3, 3) , [[5 , 5 , 5] , [5 , 5 , 5] , [ 5 , 5 , 5 ]])
		self.assertEqual(AB , A * B)
		
	def test_trace(self):
		A =  Matrix((3 , 3) , [[0 , 1 , 2 ] , [0 , 1 , 2 ] , [0 , 1 , 2 ]] )
		self.assertEqual(A.trace() , 3)
		
	def test_det(self):
		A =  Matrix((3 , 3) , [[7 , 1 , 3 ] , [1 , 1 , 4 ] , [1 , 1 , 5 ]] )
		self.assertEqual(A.determinant() , 6)
		
	def test_inverse(self):
		A =  Matrix((3 , 3) , [[1 ,1 , 2] , [1 , 1 , 1] , [2 , 1 , 0]])
		Ainv = Matrix((3 , 3) , [[1 , -2 , 1] , [-2 , 4 , -1] , [1 , -1 , 0]])
		self.assertEqual(A.inverse() , Ainv)
		
if __name__ == "__main__":
	unittest.main()
