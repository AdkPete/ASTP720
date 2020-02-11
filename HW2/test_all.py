import unittest
from matrix import *
from calc import *
import numpy as np

a = np.random.rand()
def rand_matrix(N):
	'''
	Takes in an integer N
	returns a randomly generated NxN Matrix
	'''
	A =  Matrix((N, N))
	for i in range(N):
		for j in range(N):
			A.elements[i][j] = int(np.random.rand() * 15)
			
	while A.determinant() == 0: ###If our matrix is singular, we make a new random matrix
		A =  Matrix((N, N))
		for i in range(N):
			for j in range(N):
				A.elements[i][j] = int(np.random.rand() * 15)
	return A
	
def I(N):
	'''
	Takes in an integer N
	Returns an NxN identity matrix
	'''
	
	Ident = Matrix((N , N))
	for i in range(N):
		Ident.elements[i][i] = 1
	return Ident


class TestMatrix(unittest.TestCase):


	def test_get(self):
		###tests the get function
		A = Matrix((2 , 2) , [[0 , 1 ] , [2 , 3]])
		self.assertEqual(A.get(1,1) , 0)
		
	def test_add(self):
		###Tests adding two matrices together
		A = Matrix((3 , 3) , [[0 , 1 , 2 ] , [0 , 1 , 2 ] , [0 , 1 , 2 ]] )
		B = Matrix((3 , 3) , [[0 , 1 , 2 ] , [0 , 1 , 2 ] , [0 , 1 , 2 ]] )
		Answer = Matrix((3 , 3) , [[0 , 2 , 4] , [0 , 2 , 4] , [0 , 2 , 4]])
		C = A + B
		self.assertEqual(C , Answer)
		
	def test_transpose(self):
		###Tests the transpose of a matrix
		
		A = Matrix( (3 , 3) , [[0 , 0 , 0 ] , [1 , 1 , 1 ] , [2 , 2 , 2 ]])
		AT = Matrix( (3 , 3) , [[0 , 1 , 2] , [0 , 1 , 2] , [0 , 1 , 2]])
		self.assertEqual(AT , A.transpose())
		
	def test_mult(self):
	
		###Tests multiplying two matrices together and multiplying a matrix by a scalar
		
		A =  Matrix((3 , 3) , [[0 , 1 , 2 ] , [0 , 1 , 2 ] , [0 , 1 , 2 ]] )
		B =  Matrix( (3 , 3) , [[0 , 0 , 0 ] , [1 , 1 , 1 ] , [2 , 2 , 2 ]])
		AB = Matrix( (3, 3) , [[5 , 5 , 5] , [5 , 5 , 5] , [ 5 , 5 , 5 ]])
		C = Matrix((3 , 3) , [[0 , 2 , 4 ] , [0 , 2 , 4 ] , [0 , 2 , 4 ]])
		D = A * 2
		self.assertEqual(D , C)
		self.assertEqual(AB , A * B)
		
		
	def test_trace(self):
		###Checks the trace of a matrix
		A =  Matrix((3 , 3) , [[0 , 1 , 2 ] , [0 , 1 , 2 ] , [0 , 1 , 2 ]] )
		self.assertEqual(A.trace() , 3)
		
	def test_det(self):
		###tests the determinant of a matrix
		A =  Matrix((3 , 3) , [[7 , 1 , 3 ] , [1 , 1 , 4 ] , [1 , 1 , 5 ]] )
		self.assertEqual(A.determinant() , 6)
		
	def test_inverse(self):
		### tests the inverse of a 3x3. Another inverse test is located below
		A =  Matrix((3 , 3) , [[1 ,1 , 2] , [1 , 1 , 1] , [2 , 1 , 0]])
		Ainv = Matrix((3 , 3) , [[1 , -2 , 1] , [-2 , 4 , -1] , [1 , -1 , 0]])
		self.assertEqual(A.inverse() , Ainv)
		
	def test_row_swap(self):
		###	swaps two rows
		A =  Matrix((3 , 3) , [[1 ,1 , 2] , [1 , 1 , 1] , [2 , 1 , 0]])
		A.swap_rows(0 , 1)
		C = Matrix((3 , 3) , [[1 , 1 , 1] , [1 ,1 , 2] , [2 , 1 , 0]])
		self.assertEqual(C , A)
		
	def test_LU(self):
		
		##tests lu factorizatipn
		
		A = Matrix((3 , 3) , [[1 , 2 , 3] , [ 3 , 2 , 3] , [4 , 5 , 6]])
		L , U = A.LU()
		CL = Matrix((3 , 3) , [[1 , 0 , 0 ] , [3 , 1 , 0] , [4 , .75 , 1]])
		CU = Matrix((3 , 3) , [[1 , 2 , 3] , [0 , -4 , -6] , [0 , 0 , -1.5]])
		self.assertEqual(L , CL)
		self.assertEqual(U , CU)
		
	def test_solve(self):
	
		###tests solving a system of equations with LU factorization
		A = Matrix((3 , 3) , [[3 , 3 , 3] , [-3 , 3 , 3] , [-3 , -3 , 3]])
		b = Matrix((3 , 1) , [ [9] , [9] , [9] ])
		rx = Matrix((3 , 1) , [ [ 0 ] , [ 0  ], [3 ] ])
		x = solve_eq(A , b)
		self.assertEqual(x , rx)
		
class TestCalc(unittest.TestCase):

	#runs a number of test cases

	def test_centered_diff(self):
		self.assertAlmostEqual( 5 , centered_diff(test_func_1 ,1e-4)(5))
		
	def test_rectangle(self):
		###note the low accuracy restriction due to the inefficient nature of this method
		self.assertEqual( round(right_rect_integrate(test_func_1 , 0 , 5 , 1e-5) , 3) , 112.5)
		
	def test_trap(self):
		###tests the trap rule method
		
		self.assertAlmostEqual( trap_rule(test_func_1 , 0 , 5 , 1e-4) , 112.5)
		self.assertAlmostEqual( trap_rule(test_func_2 , 0 , np.pi * 2 , 1e-4) , 0)
		
	def test_midpoint(self):
		###tests the midpoint rule methods
		self.assertAlmostEqual( midpoint_rule(test_func_1 , 0 , 5 , 1e-4) , 112.5)
		self.assertAlmostEqual( midpoint_rule(test_func_2 , 0 , np.pi * 2 , 1e-4) , 0)
		
	def test_simpson(self):
		###tests the simspons rule methods
		self.assertAlmostEqual( simpson(test_func_1 , 0 , 5 , 1e-4) , 112.5)
		self.assertAlmostEqual( simpson(test_func_2 , 0 , np.pi * 2 , 1e-4) , 0)

N = 6
A = rand_matrix(N)
b = Matrix((1 , N))

class TestRand(unittest.TestCase):

	def test_inverse(self):
	
		###Tests that A * Ainverse = I for a random matrix A
		self.assertEqual((A.inverse() * A * 5).round_all(7) , I(N) * 5)
		

		
unittest.main()
