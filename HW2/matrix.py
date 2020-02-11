import unittest
import numpy as np
import copy

class Matrix:
	'''
	A class used to perform a variety of matrix computations
	elements contains all of the entries in our matrix
	size provides the dimensions of the matrix
	'''
	
	def __init__(self , size , elements = []):
		'''
		Creates a Matrix object
		'''
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
		The one indexing is occasionally useful since many matrix equations are one indexed.
		'''
		return self.elements[i - 1][j - 1]
		
	def __eq__(self , other):
		
		'''
		We definte two matrices to be equal if every element in the two matrices is equal
		'''
		
		if self.size != other.size:
			return False
		for row in range(len(self.elements)):
			for col in range(len(self.elements[row])):
				if self.elements[row][col] != other.elements[row][col]:
					return False
		return True
		
	def print(self):
		'''
		Prints out a representation of a matrix
		Prints all elements
		'''
		
		for row in range(len(self.elements)):
			row_out = ""
			for col in range(len(self.elements[row])):
				row_out += str(self.elements[row][col]) + "\t\t"
			print (row_out)
	
	def swap_rows(self, i , j):
		'''
		This function will swap the locations of two rows, whose indices are i and j
		Note that these indices are zero indexed.
		'''
		
		row_i = self.elements[i]
		row_j = self.elements[j]
		self.elements[i] = row_j
		self.elements[j] = row_i
		
	def round_all(self , n):
		'''
		takes in an integer n
		rounds every element to n decimal places
		useful for testing purposes
		'''
		New = Matrix(self.size)
		for i in range(len(self.elements)):
			for j in range(len(self.elements[i])):
				New.elements[i][j] = round(self.elements[i][j] , n)
		return New
		
	def __add__(self , other):
	
		'''
		Adds two matrices together
		Only works if given two Matrix objects
		returns a new Matrix object
		'''
		
		New = Matrix(self.size)
		if self.size != other.size:
			 raise SystemExit('Attempted to add matrices that are not the same size')
		for row in range(len(self.elements)):
			for col in range(len(self.elements[row])):
				New.elements[row][col] = self.elements[row][col] + other.elements[row][col]
		return New

	def transpose(self):
		'''
		Takes the transpose of a matrix
		returns a new Matrix object
		'''
		New = Matrix((self.size[0] , self.size[1]))
		for row in range(len(self.elements)):
			for col in range(len(self.elements[row])):
				New.elements[col][row] = self.elements[row][col]
				
		return New
		
	def __mul__(self , other):
		
		'''
		If other is another matrix, we perform matrix multiplication
		If other is a scalar, then we multiply each element by said scalar
		Returns a new Matrix object
		'''
		
		if type(other) == type(1) or type(other) == type(5.0):
			N = Matrix(self.size)
			for i in range(len(N.elements)):
				for j in range(len(N.elements[i])):
					N.elements[i][j] = self.elements[i][j] * other
			return N
		N = Matrix((self.size[0] , other.size[1]))
		
		for row in range(len(N.elements)):
			for col in range(len(N.elements[row])):
				V = 0
				for i in range(len(self.elements)):
					V += self.elements[row][i] * other.elements[i][col]
				
				N.elements[row][col] = V
		
		return N
		
	def trace(self):
	
		'''
		Will return the trace of a function
		returns the value of the trace
		'''
		
		TR = 0
		if self.size[0] != self.size[1]:
			raise SystemExit('Cannot calculate the trace of a non-square matrix')
		for i in range(self.size[0]):
			TR += self.elements[i][i]
		
		return TR
		
		
	def residual(self , i , j):
		'''
		This function determines the residual matrix given a particular element specified by i and j
		i and j are zero indexed here
		This function primarily exists because it is useful for taking determinants and inverses
		returns a new Matrix object
		'''
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
		'''
		Calculates the determinant of a matrix
		Operates recursively using residual matrices
		returns the value of the determinant
		'''
		
		if self.size == (2 , 2):
			return self.elements[0][0] * self.elements[1][1] - self.elements[0][1] * self.elements[1][0]
		###Now we handle larger matrices
		det = 0
		for i in range(len(self.elements)):
			Res = self.residual(0 , i)
			det += (-1) ** (i + 2) * self.elements[0][i] * Res.determinant()
		
		return float(det)
		
	def inverse(self):
	
		'''
		Calculates the inverse of a matrix
		returns a new matrix object
		'''
		
			
		New = Matrix(self.size)
		D = self.determinant()
		if D == 0:
			raise SystemExit('Attempted to invert a singular matrix')
		
		if self.size == (2 , 2): ###handles inverting a 2x2
		
			New.elements[0][0] = self.elements[1][1]
			New.elements[1][1] = self.elements[0][0]
			New.elements[0][1] = -1 * self.elements[1][0]
			New.elements[1][0] = -1 * self.elements[0][1]
			return New
	
		
		for row in range(len(self.elements)):
			for col in range(len(self.elements[row])):
				Cij = (-1) ** (row + col + 2) * self.residual(row , col).determinant()
				New.elements[row][col] = Cij / D
		return New.transpose()
		
		
	def LU(self):
	
		'''
		Performs an LU decomposition on the Matrix
		Will return two new Matrix objects, L and U
		'''
		L = Matrix(self.size)
		U = Matrix(self.size)
		
		for i in range(len(L.elements)):
			for j in range(i , len(L.elements[i])):
				sum = 0
				k = 0
				while k < i:
					sum += L.elements[i][k] * U.elements[k][j]
					k += 1
				U.elements[i][j] = self.elements[i][j] - sum
				
			for j in range(i , len(L.elements[i])):
				if i == j:
					L.elements[i][j] =1
					continue
				sum = 0
				k = 0
				while k < i:
					sum += L.elements[j][k] * U.elements[k][i]
					k += 1
				L.elements[j][i] = (self.elements[j][i] - sum) / U.elements[i][i]
			
		return L , U
		
		

def solve_eq(A , b):

	'''
	Takes in two matrix objects, A and b.
	Solves a system of equations given by Ax = b
	returns a Matrix object x
	'''
	
	L , U = A.LU()
	y = [0] * len(L.elements)
	y[0] = b.elements[0][0] / L.elements[0][0]
	
	
	for i in range(1 , len(y)):
		
		k = 0
		Sum = 0
		while k < i:
			Sum += L.elements[i][k] * y[k]
			k += 1
		y[i] = (1 / L.elements[i][i]) * (b.elements[i][0] - Sum)
		
	
	x = [0] * len(y)
	x[-1] = y[-1] / U.elements[-1][-1]
	for i in range( len(y) - 2, -1 , -1):
		Sum = 0
		k = i
		
		while k < len(y):
			Sum += U.elements[i][k] * x[k]
			k += 1
			
		x[i] = (1 / U.elements[i][i]) * (y[i] - Sum)
	
	return Matrix((len(x) , 1) , [ x ])
